"""
Created on May 06, 2014

@author: william

    This script creates different libraries templates.
"""

from __future__ import division
import ast
import ConfigParser
from ConfigParser import NoOptionError
import os
from multiprocessing import Pool
import sys

from astropy import units
import h5py
import numpy as np
from magal.index.lick import LickIndex

from magal.library import LibraryModel
from magal.io.readfilterset import FilterSet
from magal.util.cosmo import zcor
from magal.photometry.syntphot import spec2filterset


def spec2filter_z(args):  # FIXME: Move this to somewhere else
    spec, ccd_filter, z, cosmo = args
    n_spec = np.shape(spec)[0]

    if z == 0:
        d_L = 3.08567758e19  # 10 parsec in cm
    else:
        d_L = cosmo.luminosity_distance(z).to('cm')
    k_cosmo = units.Lsun.to('erg/s') / (4 * np.pi * np.power(d_L, 2))

    if n_spec == 1:  # If there is only model.
        model_spec = zcor(spec[0], z)
        model_spec['flux'] = model_spec['flux'] * k_cosmo
        x = spec2filterset(ccd_filter, model_spec)
    elif n_spec == 2:
        obs_spec = zcor(spec[0], z)
        obs_spec['flux'] = obs_spec['flux'] * k_cosmo
        obs_spec['error'] = obs_spec['error'] * k_cosmo
        model_spec = zcor(spec[1], z)
        model_spec['flux'] = model_spec['flux'] * k_cosmo
        x = spec2filterset(ccd_filter, obs_spec, model_spec)

    return x


def main():
    if len(sys.argv) < 2:
        print 'Usage: %s configuration.ini' % sys.argv[0]
        return 1

    # ### ALL the configfile variables are loaded here ####
    # TODO: Put some checking/error handling here.
    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])
    libfile = os.path.expandvars(config.get('LibraryGeneral', 'library_file'))
    library_type = os.path.expandvars(config.get('LibraryGeneral', 'library_type'))
    filter_file = os.path.expandvars(config.get('LibraryGeneral', 'filter_file'))

    filter_id = config.get('LibraryGeneral', 'filterset')  # TODO: Check if it is valid.
    z_from = config.getfloat('LibraryGeneral', 'z_from')
    z_to = config.getfloat('LibraryGeneral', 'z_to')
    z_step = config.getfloat('LibraryGeneral', 'z_step')
    try:
        z_error = config.getfloat('LibraryGeneral', 'z_error')
        print 'z_error is ', z_error
    except NoOptionError:
        z_error = None
        print 'z_error is None'
    try:
        cosmology = ast.literal_eval(config.get('LibraryGeneral', 'cosmology'))
        from astropy.cosmology import FlatLambdaCDM

        cosmo = FlatLambdaCDM(H0=cosmology['H0'], Om0=cosmology['Om0'])
    except NoOptionError:
        print 'Adopting the WMAP9 cosmology:'
        from astropy.cosmology import WMAP9 as cosmo

        print cosmo.__doc__

    try:
        allow_overwrite = config.getboolean('LibraryGeneral', 'allow_overwrite')
    except NoOptionError:
        allow_overwrite = False
    if os.path.exists(libfile):
        if allow_overwrite:
            try:
                os.unlink(libfile)
            except:
                print u'Cannot DELETE file {0:s}.'.format(libfile)
                return 1
        else:
            print u'File {0:s} already exists.'.format(libfile)
            return 1


    # 1 - Two-exponential library parameters:
    if library_type == 'two_exp':
        inp_file = os.path.expandvars(config.get('LibraryParameters', 'input_file'))
        basedir = os.path.expandvars(config.get('LibraryParameters', 'bases_dir'))
        basefile = os.path.expandvars(config.get('LibraryParameters', 'base_file'))
        fraction_type = config.get('LibraryParameters', 'fraction_type')  # light or mass
        try:
            lambda_norm = config.getfloat('LibraryParameters', 'lambda_norm')  # mass normalization wl
        except NoOptionError:
            lambda_norm = 4020
        LibModel = LibraryModel(library_type, inp_file, basedir, basefile, fraction_type, lambda_norm)  # Init
    # 2 - Two-exponential library w/ random parameters
    elif library_type == 'two_exp_random':
        inp_file = os.path.expandvars(config.get('LibraryParameters', 'input_file'))
        basedir = os.path.expandvars(config.get('LibraryParameters', 'bases_dir'))
        basefile = os.path.expandvars(config.get('LibraryParameters', 'base_file'))
        fraction_type = config.get('LibraryParameters', 'fraction_type')  # light or mass
        try:
            lambda_norm = config.getfloat('LibraryParameters', 'lambda_norm')  # mass normalization wl
        except NoOptionError:
            lambda_norm = 4020
        LibModel = LibraryModel(library_type, inp_file, basedir, basefile, fraction_type, lambda_norm)  # Init
    # 3 - STARLIGHT library parameters:
    elif library_type == 'starlight_sdss':
        library_type = 'starlight_sdss'
        inp_file = os.path.expandvars(config.get('LibraryParameters', 'input_file'))
        tables_dir = os.path.expandvars(config.get('LibraryParameters', 'tables_dir'))
        input_dir = os.path.expandvars(config.get('LibraryParameters', 'input_dir'))
        output_dir = os.path.expandvars(config.get('LibraryParameters', 'output_dir'))
        LibModel = LibraryModel(library_type, inp_file, tables_dir, input_dir, output_dir, cosmo)

    #### -- ####

    # 0 - Library Parameters
    z_range = np.arange(z_from, z_to, z_step)
    Nz = len(z_range)

    # 1 - Init/Read Files
    ## 1.1 - Load filter information
    f = FilterSet(filter_file)
    ## 1.2 - New Library
    db = h5py.File(libfile, 'w')

    # Write full .ini file on model outputfile.
    aux_ini = open(sys.argv[1], 'r').read()
    ds = db.create_dataset('/ini_file', shape=(1,), dtype=h5py.new_vlen(str))
    ds[:] = aux_ini

    ### 1.2.1 - Tables group
    db.create_group('/tables/')
    ### 1.2.2 - Filtersystem group.
    db.create_group('/%s/' % filter_id)
    ### 1.2.3 - CCD groups
    ccd_vec = {}
    db_vec = {}
    for ccd in f.filtersets[filter_id]:
        f.load(filter_id, ccd)  # Load CCD filters.
        filter_wls = f.filter_wls  # Get the wls
        Nl = len(filter_wls)
        db.create_group('/%s/%s/' % (filter_id, ccd))  # Create group to the CCD.

        ccd_vec.update({ccd: f.filterset})

        db.create_dataset(name='/%s/%s/filtercurves' % (filter_id, ccd),
                          data=f.filterset)  #np.array(f.filterset, dtype=([('ID_filter', '|S32'), ('wl', '<f4'), ('transm', '<f4')])))
        db.create_dataset(name='/%s/%s/filterset' % (filter_id, ccd), data=filter_wls)
        chunk_size = 1 * LibModel.lib_size * Nl  # chunk_size = 1 (z) * lib_size * n_filters
        # chunk_z = int(chunk_size / (64*1024) ) + 1 # make chunk_size of about 64kb and > 1
        chunk_z = 1  # print (Nz,lib_size, Nl), (chunk_z,lib_size,Nl)
        db_m = db.create_dataset('/%s/%s/%s' % (filter_id, ccd, 'library'), shape=(Nz, LibModel.lib_size, Nl),
                                 chunks=(chunk_z, LibModel.lib_size, Nl),
                                 dtype=np.dtype([('m_ab', '<f4'), ('e_ab', '<f4')]))
        db_vec.update({ccd: db_m})
    ### 1.2.3 - redshift
    db.create_dataset(name='/tables/z', data=z_range)
    if z_error:
        z_spec = np.array(
            [np.random.normal(loc=z_range[i_z], size=LibModel.lib_size, scale=(1 + z_range[i_z]) * z_error) for i_z in
             range(Nz)],
            dtype=np.float)
        z_spec[z_spec < 0.001] = 0.001  # Avoid zeroes
        db.create_dataset(name='/tables/z_spec', data=z_spec)
    else:
        z_spec = None
    ### 1.2.4 - properties
    prop = db.create_dataset(name='/tables/properties', data=LibModel.input_data)

    ### 1.2.5 - indexes
    try:
        conf_indexes = ast.literal_eval(config.get('LibraryGeneral', 'indexes'))
        conf_indexes.sort()
        N_indexes = len(conf_indexes)
        dt = np.dtype([('index', '<f4'), ('error', '<f4')])
        indexes = np.empty(shape=(LibModel.lib_size, N_indexes), dtype=dt)
        Index = LickIndex()
        calc_index = True
    except NoOptionError:
        print 'caca'
        calc_index = False
        pass



    # 2 - RUN!
    for i_model in range(LibModel.lib_size):
        if i_model % 100 == 0:
            print 'Running i_model: %i' % i_model
        for ccd in f.filtersets[filter_id]:
            spec = LibModel.get_model_spectrum(i_model)
            if calc_index:
                aux = Index.lick_index(conf_indexes, spec[1]['wl'], spec[1]['flux'])
                for i_index in range(len(conf_indexes)):
                    indexes[i_model, i_index]['index'] = np.float32(aux[np.sort(conf_indexes)[i_index]][0])
                    indexes[i_model, i_index]['error'] = np.float32(aux[np.sort(conf_indexes)[i_index]][1])
                calc_index = False
            if z_error:
                aux_z = z_spec[:, i_model]
            else:
                aux_z = z_range
            args = [(spec, ccd_vec[ccd], aux_z[i_z], cosmo) for i_z in range(Nz)]
            pool = Pool()
            result = pool.map(spec2filter_z, args)
            pool.close()
            pool.join()
            db_vec[ccd][:, i_model, :] = result

    # caca = db.create_dataset(name='/tables/spectral_indexes', data=indexes)
    # caca.attrs.create(name='index_names', data=conf_indexes)

    print 'Finished calculating magnintudes.'

    db.close()