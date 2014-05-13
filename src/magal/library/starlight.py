import gzip
import re
from astropy import units
import bz2

from magal.core.exceptions import MAGALException
from magal.util.matchs import matchobjs

import numpy as np
import h5py
import atpy
import pystarlight.io

__author__ = 'william'


class StarlightSDSS(object):
    def __init__(self, type, filename, tables_dir, input_dir, output_dir, data_release='DR7', sample='926246',
                 base='BS'):
        self.type = type
        self._input_dir = input_dir
        self._output_dir = output_dir

        self._sl_files = np.loadtxt(filename, dtype=np.dtype([('input', 'S32'), ('output', 'S32')]))
        self.lib_size = len(self._sl_files)
        self._el_file = '%s/sample.F%s.%s.f.lines.dat.%s.hdf5' % (tables_dir, data_release, sample, base)
        self._syn01_file = '%s/sample.F%s.%s.f.Starlight.SYN01.tab.%s.hdf5' % (tables_dir, data_release, sample, base)
        self._syn02_file = '%s/sample.F%s.%s.f.Starlight.SYN02.tab.%s.hdf5' % (tables_dir, data_release, sample, base)
        self._syn03_file = '%s/sample.F%s.%s.f.Starlight.SYN03.tab.%s.hdf5' % (tables_dir, data_release, sample, base)
        self._syn04_file = '%s/sample.F%s.%s.f.Starlight.SYN04.tab.%s.hdf5' % (tables_dir, data_release, sample, base)


    @property
    def input_data(self):
        return self._get_data()

    def _get_data(self):
        ## SYN [1-4] properties
        syn = {}
        syn['syn01'] = h5py.File(self._syn01_file, 'r').get('/SYN01')
        syn['syn02'] = h5py.File(self._syn02_file, 'r').get('/SYN02')
        syn['syn03'] = h5py.File(self._syn03_file, 'r').get('/SYN03')
        syn['syn04'] = h5py.File(self._syn04_file, 'r').get('/SYN04')
        syn_prop = {'syn01': ['A_V', 'v0', 'vd', 'SN_w', 'SN_n'],
                    'syn02': ['at_flux', 'at_mass', 'aZ_flux', 'aZ_mass', 'am_flux', 'am_mass'],
                    'syn03': ['M2L_r'],
                    'syn04': ['Mcor_fib', 'Mcor_gal', 'z']}  # 'DL_Mpc', 'Mini_fib' comes from SYN03

        ## Emission lines
        ### The emission lines is a bit tricky. First I separate the items from the file I want, then I calculate some
        ### emission line indexes that would be useful. I don't imagine how to do this without hardcoding.
        re_elines = re.compile('^EW|^F')
        tb_el = h5py.File(self._el_file, 'r').get('emissionLines')
        el_names = []
        for n_ in tb_el.dtype.names:
            if re_elines.match(n_):
                el_names.append(n_)

        ## Emission line ratios
        ### This, again, could be un-hardcoded.
        el_ratios = {'N2Ha': [['F_6584'], ['F_6563']],
                     'O3Hb': [['F_5007'], ['F_4861']],
                     'HaHb': [['F_6563'], ['F_4861']],
                     'O2Hb': [['F_3727'], ['F_4861']],
                     'O3N2': [['F_5007'], ['F_6584']],
                     'S2Ha': [['F_6717', 'F_6731'], ['F_6563']]}

        ## Construct the final array:
        syn_names = reduce(lambda x, y: x + y, syn_prop.values())
        id_list = matchobjs(self._sl_files['output'], syn['syn01']['id'])  # Match inputfile objectlist to Starlight tb.
        n_galaxies = len(id_list)
        dt = [(name, np.float) for name in syn_names + el_names + el_ratios.keys()]
        dt.insert(0, ('id', 'S32'))  # file-id, from SYN01 table
        data = np.empty(n_galaxies, dtype=dt)

        ## Populate the final array
        ### SYN [1-4]
        data['id'] = syn['syn01']['id'][id_list]
        for s_file in syn_prop.keys():
            for p_ in syn_prop[s_file]:
                data[p_] = syn[s_file][p_][id_list]
        ### Emission Lines
        for p_ in el_names:
            mask = (tb_el[p_][id_list] < -998)
            data[p_] = np.log10(tb_el[p_][id_list])
            data[p_][mask] = -999
        ### Emission line ratios:
        for p_ in el_ratios.keys():
            a = np.zeros(n_galaxies)
            b = np.zeros(n_galaxies)
            mask = np.zeros(n_galaxies, dtype=bool)
            for el in el_ratios[p_][0]:
                a += tb_el[el][id_list]
                mask += (tb_el[el][id_list] <= -998)
            for el in el_ratios[p_][1]:
                b += tb_el[el][id_list]
                mask += (tb_el[el][id_list] <= -998)
            data[p_] = np.log10(a / b)
            data[p_][mask] = -999

        ## Get aux_units from input tables to correct the units of the model/observed STARLIGHT spectra.
        self._aux_units = 4 * np.pi * np.power(syn['syn04']['DL_Mpc'][id_list] * units.Mpc.to('cm'), 2) / (
            np.power(10, syn['syn04']['Mini_fib'][id_list]) * units.Lsun.to('erg/s'))

        ## Close all files
        tb_el.file.close()
        for f_ in syn.keys():
            syn[f_].file.close()

        ## Mask invalid columns w/ nan.
        mask = [tuple(np.append(False, np.less(list(p_)[1:], -998))) for p_ in data]
        return np.ma.masked_array(data, mask=mask, fill_value=np.nan).filled()

    def get_model_spectrum(self, i_spec):
        """
        Return the spectrum of a given model on the models list.

        Parameters
        ----------
        i_spec : int
            i_spec number which we want to take the spectrum.

        Returns
        -------
        spec : array
            Two arrays with ``lambda`` and ``flux`` to the observed and model ``i_spec`` spectrum.
        """
        # Read starlight output
        model_file = '%s/%s' % (self._output_dir, self._sl_files['output'][i_spec])
        tm = atpy.TableSet(model_file, type='starlightv4')
        model_spec = np.copy(tm.spectra.data.view(dtype=np.dtype(
            [('wl', tm.spectra.l_obs.dtype), ('f_obs', '<f8'), ('flux', tm.spectra.f_syn.dtype), ('f_wei', '<f8'),
             ('Best_f_SSP', '<f8')])))
        #     and input..
        obs_file = '%s/%s' % (self._input_dir, self._sl_files['input'][i_spec])
        if obs_file.endswith('.gz'):
            f = gzip.GzipFile(obs_file)
        elif obs_file.endswith('.bz2'):
            f = bz2.BZ2File(obs_file)
        else:
            f = open(obs_file)
        obs_spec = np.loadtxt(f, dtype=np.dtype([('wl', '<f8'), ('flux', '<f8'), ('error', '<f8'), ('flag', '<i8')]))

        # If aux_units not defined, will throw a NameError exception. This is expected because you did not ran
        # input_data() before doing this...
        model_spec['flux'] = model_spec['flux'] * tm.keywords['fobs_norm'] * 1e-17 * self._aux_units[i_spec]
        obs_spec['flux'] = obs_spec['flux'] * 1e-17 * self._aux_units[i_spec]
        obs_spec['error'] = obs_spec['error'] * 1e-17 * self._aux_units[i_spec]

        return obs_spec, model_spec

    def _check_input(self):
        if self.type != 'starlight_sdss':
            print('Error creating StarlightSDSS LibraryModel object.')
            raise MAGALException()
        pass