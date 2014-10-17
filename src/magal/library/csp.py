__author__ = 'william'

import os

import atpy
from pystarlight.util.redenninglaws import Charlot_RedLaw
import numpy as np

from ..core.exceptions import MAGALException
from ..util.stellarpop import n_component


class TwoExponential(object):
    """
    A simple library w/ two exponential decaying components.
    """

    def __init__(self, library_type, input_file, bases_dir, base_file, fraction_type, lambda_norm=4020.):
        """
        A simple library w/ two exponential decaying components.

        Parameters
        ----------
        library_type : string
            Will raise an error if it is different from ``two_exp``.

        input_file : string
            Base models input filename. The file must be an ascii file with this 7 ``float`` the columns:
            ``t0_young, tau_young, t0_old, tau_old, frac_young, tau_v, Z``

        bases_dir : string
            Directory where the models files are.

        base_file : string
            Starlight-like base filename.

        fraction_type : string
            Type of weighting the young component ``frac_young``. Can be ``light`` or ``mass`` weighting.

        lambda_norm : float, optional
            If ``fraction_type`` is set to be on ``mass``, it will use this as mass-normalization lambda.
        """

        # 0 - Check if all the parameters are okay.
        self.type = library_type
        self.input_file = input_file
        self.bases_dir = bases_dir
        self.base_file = base_file
        self.fraction_type = fraction_type  #FIXME: Add check
        self.lambda_norm = lambda_norm
        self._check_input()

        # 1 - Load data...
        # 1.1 - Input File w/ library parameters
        dt = np.dtype([('t0_young', np.float), ('tau_young', np.float), ('t0_old', np.float), ('tau_old', np.float),
                       ('frac_young', np.float), ('tau_v', np.float), ('Z', np.float)])
        self.input_data = np.loadtxt(self.input_file, dtype=dt)
        self.lib_size = len(self.input_data)

        # 1.2 - STARLIGHT base file + dir.
        self.bt = atpy.Table(self.base_file, self.bases_dir, read_basedir=True, type='starlightv4_base')
        self.ages = np.unique(self.bt['age_base'])

        self._wl = self.bt[0]['l_ssp']  # FIXME: Put some error handling when we have ssps w/ different wl coverages.

        if self.fraction_type == 'mass':
            self._i_norm = np.argmin((self._wl - self.lambda_norm) ** 2)

        # 1.3 - Charlot & Fall 2000 reddening law. FIXME: Abstract this to use other extinction laws.
        self.tau_l_Y, self.tau_l_O = Charlot_RedLaw(self._wl,
                                                    mu=0.3)  # FIXME: Put some error handling when we have ssps w/ different wl coverages.


    def get_model_spectrum(self, i_model):
        """
        Return the spectrum of a given model on the models list.

        Parameters
        ----------
        i_model : int
            Model number which we want to take the spectrum.

        Returns
        -------
        spec : array
            An array with ``lambda`` and ``flux`` to the model ``i_model`` spectrum.
        """

        # 0 - Init some vars and flags
        spec = np.zeros(shape=np.shape(self.bt[0]['f_ssp']), dtype=np.dtype([('wl', np.float), ('flux', np.float)]))
        spec['wl'] = self._wl
        flag_Z = (self.bt['Z_base'] == self.bt['Z_base'][np.argmin(
            (self.input_data[i_model]['Z'] - self.bt['Z_base']) ** 2)])  # Get a mask to the desired metallicity

        # 1 - Init the n_component model and the output spectra
        model = n_component(self.ages)

        # 2 - Add the components to our model
        # 2.1 - Young
        model.add_exp(self.input_data[i_model]['t0_young'], self.input_data[i_model]['tau_young'], 1)
        light_fraction = self.input_data[i_model]['frac_young']
        # 2.2 - Old
        if self.input_data[i_model]['frac_young'] < 1:
            model.add_exp(self.input_data[i_model]['t0_old'], self.input_data[i_model]['tau_old'], 1)
            light_fraction = [light_fraction, 1 - self.input_data[i_model]['frac_young']]
        # sfh = model.get_sfh() #TODO: Why?

        spec_components = []
        l_norm = []
        for sfh_curve in model.individual_sfh_curves:
            aux_spec = np.zeros_like(spec['flux'])
            for i in range(len(sfh_curve)):
                if sfh_curve[i] > 0:
                    if model.ages_start[i] <= 1e7:  # We divide the evaluation of the extinction law coefficients to
                        # match Charlot&Fall 2000 method.
                        aux_spec += self.bt[flag_Z][i]['f_ssp'] * sfh_curve[i] * np.exp(
                            -self.input_data[i_model]['tau_v'] * self.tau_l_Y) / self.bt[i]['Mstars']
                    else:
                        aux_spec += self.bt[flag_Z][i]['f_ssp'] * sfh_curve[i] * np.exp(
                            -self.input_data[i_model]['tau_v'] * self.tau_l_O) / self.bt[i]['Mstars']
            l_norm.append(aux_spec[self._i_norm])
            spec_components.append(aux_spec)

        if self.fraction_type == 'mass':  # If frac_young is set to be as mass fraction, then correct it.
            mass_fraction = np.array(light_fraction) / np.array(l_norm)
            mass_fraction /= np.sum(mass_fraction)
            for i in range(len(model.individual_sfh_curves)):
                spec['flux'] += spec_components[i] * mass_fraction[i]
        elif self.fraction_type == 'light':
            spec['flux'] = aux_spec

        return spec[np.newaxis,:]  # For CSPs, we return only the model spectrum

    def _check_input(self):
        if self.type != 'two_exp':
            print('Error creating TwoExponential LibraryModel object.')
            raise MAGALException()

        for f_ in [self.input_file, self.base_file]:
            if os.path.isfile(f_):
                pass
            else:
                print('File %s not found or is not a file.' % f_)
                raise MAGALException('File %s not found or is not a file.' % f_)
        if os.path.isdir(self.bases_dir):
            pass
        else:
            print('Directory %s not found or is not a directory.' % self.bases_dir)
            raise MAGALException('File %s not found or is not a file.' % f_)

