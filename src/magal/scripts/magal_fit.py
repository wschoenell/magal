#!/usr/bin/env python
# encoding: utf-8
'''
MAGAL -- Magnitudes Analyzer fitting program.

@author:     william
        
@license:    GPLv3

@contact:    william@iaa.es
'''
import hashlib
import logging
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import multiprocessing

import sys
import h5py
import numpy as np
import time
from magal.core.config import MagalConfig

from magal.core.log import logger
from magal.io.readlibrary import Library
from magal.fit.stats import percentiles
from magal.io.readinput import Input
from magal.util.matchs import get_zslice
from magal.fit.magal_aux import chi2_parameters, chi2_wrapper

DEBUG = 0
TESTRUN = 0
PROFILE = 0


def md5_for_file(f, block_size=2 ** 20):
    fp = open(f, 'rb')
    md5 = hashlib.md5()
    while True:
        data = fp.read(block_size)
        if not data:
            break
        md5.update(data)
    fp.close()
    return md5.hexdigest()


def main(argv=None):  # IGNORE:C0111
    log = logger(__name__)
    logging.root.setLevel(logging.DEBUG)
    log.setLevel(logging.DEBUG)

    ini_file = sys.argv[1]
    config = MagalConfig(ini_file, 'fit')
    if config.ret > 0:
        return config.ret

    t0 = time.time()

    # - Init outputfile
    try:
        f = h5py.File(config.output_file, mode='w-')
    except IOError:
        print(u'File {0:s} already exists or could not be created.'.format(config.output_file))
        return 1

    # 2.2 - Define some auxiliar data...
    f.attrs.create('input_file', config.input_file)
    f.attrs.create('input_md5', md5_for_file(config.input_file))  # Store Inputfile MD5SUM
    f.attrs.create('library_file', config.library_file)
    f.attrs.create('library_md5', md5_for_file(config.library_file))  # Store Library MD5SUM

    # Write full .ini file on model outputfile.
    aux_ini = open(ini_file).read()
    ds = f.create_dataset('/ini_file', shape=(1,), dtype=h5py.new_vlen(str))
    ds[:] = aux_ini
    ###

    # 1 - Load files
    # 1.1 - Inputfile
    if config.is_simulation:  # If this is a simulation, a Library will be the inputfile.
        inp = Library(config.input_file)
    else:
        inp = Input(config.input_file)

    if config.ccd is None:
        ccds = inp.ccds
    else:
        ccds = [config.ccd]

    # Run fit over all CCDs...
    for ccd in ccds:
        log.debug('Running magal_fit on %s/%s' % (config.filter_sys, ccd))
        inp.get_filtersys(config.filter_sys, ccd)

        # 1.2 - Library
        lib = Library(config.library_file)
        lib.get_filtersys(config.filter_sys, ccd)

        # 1.3 - Get objectlist
        if config.is_simulation:
            o_list = get_zslice(inp, config.obj_z)  # Need to get selected redshift from the library
        else:
            o_list = inp.data

        if config.Nobj_max > 0:
            o_list = o_list[:config.Nobj_max]
            log.info('Running magal fit to ONLY %i objects.' % config.Nobj_max)
        else:
            log.info('Running magal fit to ALL objects.')

        N_obj = len(o_list)

        if config.is_simulation:
            if config.mass_field is not None:
                config.mass_field = inp.properties[config.mass_field][:N_obj]  # Multiply by the mass_field.

        # 2 - Output file
        # 2.3 - Data matrices
        # 2.3.1 - Shape is defined by (N_obj, N_z, N_lib)
        if config.Nz or (config.is_simulation and config.obj_z):
            N_z, N_tpl = (1, lib.library.shape[1])  # If the redshift comes from inputfile, N_z = 1.
            log.debug('Fitting only over ONE redshift!')
        if config.is_simulation and config.obj_z:
            inp_z = [config.obj_z]
        elif config.Nz:
            inp_z = inp.z
        else:
            N_z, N_tpl = lib.library.shape[:2]
            inp_z = lib.z
        aux_shape = (N_obj, N_z, N_tpl)

        log.debug('Model shape is (%i, %i, %i) - (N_obj, N_z, N_lib)' % aux_shape)

        # 2.3.2 - Create datasets
        n_ds = f.create_dataset('%s/n' % lib.path, shape=aux_shape, dtype=np.int, compression='gzip',
                                compression_opts=4)  # Number of used pixels on \chi2 calc
        s_ds = f.create_dataset('%s/s' % lib.path, shape=aux_shape, compression='gzip',
                                compression_opts=4)  # Scaling-factor
        chi2_ds = f.create_dataset('%s/chi2' % lib.path, aux_shape, compression='gzip', compression_opts=4)  # \Chi2
        mod_stats_grp = f.create_group(
            '%s/statistics/model' % lib.path)  # Model likelihood statistics (i.e. average, percentiles, etc..)
        lib_stats_grp = f.create_group('%s/statistics/library' % lib.path)  # Library likelihood statistics

        # 3 - Fit
        if config.filters_include is not None:
            filterset_mask = np.sum([k == lib.filterset['ID_filter'] for k in config.filters_include], axis=0, dtype=np.bool)
            if filterset_mask.sum() != len(config.filters_include):
                log.error(
                    'Unable to include all filters of filter_include Config Option.\n\tfilters are: %s\n' +
                    '\tfilters_include are: %s' % (lib.filterset['ID_filter'], config.filters_include))
                return 1
        elif config.filters_exclude is not None:
            filterset_mask = ~np.sum([k == lib.filterset['ID_filter'] for k in config.filters_exclude], axis=0, dtype=np.bool)
            if (~filterset_mask).sum() != len(config.filters_exclude):
                log.error(  #FIXME: TypeError: not all arguments converted during string formatting.
                    'Unable to include all filters of filter_exclude Config Option.\n\tfilters are: %s\n' +
                    '\tfilters_exclude are: %s' %
                    (', '.join(lib.filterset['ID_filter']), ', '.join(config.filters_exclude)))
                return 1
        else:
            filterset_mask = np.ones(len(lib.filterset), dtype=np.bool)

        # 3.3 - Eval \chi^2
        if config.is_simulation:
            params = chi2_parameters(config.Nz, config.is_simulation, o_list, N_obj, N_z, N_tpl, inp_z, filterset_mask,
                                     lib, config.mass_field, config.zp_error)
        else:
            params = chi2_parameters(config.Nz, config.is_simulation, o_list, N_obj, N_z, N_tpl, inp_z, filterset_mask,
                                     lib, None, config.zp_error)
        n_cpu = multiprocessing.cpu_count()
        if n_cpu > N_obj:
            map_function = map
            log.info('There is less objects than CPUs. Not using multiprocessing...')
        elif config.allow_multiprocessing:
            pool = multiprocessing.Pool()
            map_function = pool.imap
            log.debug('Started pool with %i processors.' % n_cpu)
        else:
            map_function = map
            log.info('Multiprocessing set to FALSE.')

        for result in map_function(chi2_wrapper, params):
            i_obj, aux_n, aux_s, aux_chi2 = result
            n_ds[i_obj], s_ds[i_obj], chi2_ds[i_obj] = aux_n, aux_s, aux_chi2

        # 4 - From \chi^2, eval likelihood
        # 4.1 - Likelihood
        log.debug('Calculating Likelihood. t = %3.2f' % (time.time() - t0))
        if config.cooking_factor:
            l = np.exp(-0.5 * config.cooking_factor * (np.subtract(chi2_ds, np.min(chi2_ds, axis=2).reshape((N_obj, N_z, 1)))))
        else:
            l = np.exp(-0.5 * (np.subtract(chi2_ds, np.min(chi2_ds, axis=2).reshape((N_obj, N_z, 1)))))

        # 4.2 - Normalization
        l /= np.sum(l, axis=2).reshape((N_obj, N_z, 1))  #FIXME: Double sum here!
        likelihood_ds = f.create_dataset('%s/likelihood' % lib.path, data=l, compression='gzip', compression_opts=4)

        l_T = np.sum(l, axis=1)  # Template likelihood
        norm = np.sum(l_T, axis=1)  # Template likelihood normalization
        l_T = np.array([l_T[i_obj] / norm[i_obj] for i_obj in range(N_obj)])

        f.create_dataset('%s/likelihood_template' % lib.path, data=l_T, compression='gzip', compression_opts=4)
        if N_z > 1:  # If we fit over the redshift space, calulate the redshift-compressed likelihood
            l_z = np.sum(l, axis=2)
            norm = np.sum(l_z, axis=1)
            l_z = np.array([l_z[i_obj] / norm[i_obj] for i_obj in range(N_obj)])
            lz_ds = f.create_dataset('%s/likelihood_redshift' % lib.path, data=np.array(l_z), compression='gzip',
                                     compression_opts=4)  # @UnusedVariable

        # 5 - From likelihood, eval statistics...
        # 5.0 - Some definitions...
        perc = [2.5, 16, 50, 84, 97.5]  # Percentiles which we calculate.
        perc_dt = np.dtype([('%s' % p, np.float) for p in perc])
        n_perc_dt = np.dtype([('%s' % p, np.int16) for p in perc])

        # 5.1 - Model statistics
        # 5.1.0 - Number of galaxies on each percentile.
        n_perc = mod_stats_grp.create_dataset('n_percentiles', shape=(N_obj,), dtype=n_perc_dt)
        for p_ in perc:
            n_perc[str(p_)] = np.greater_equal(np.cumsum(np.sort(l_T), axis=-1), p_/100.).sum(axis=-1)

        # 5.1.1 - Scaling-factor "s".
        log.debug('Calculating scaling-factor statistics. t = %3.2f' % (time.time() - t0))
        p_grp = mod_stats_grp.create_group('s_Mass')  # Scale prop group
        s_Mass = np.divide(s_ds.value.copy(), -2.5)

        # 5.1.2 - Best match
        log.debug('Calculating Best Matches. t = %3.2f' % (time.time() - t0))
        bmx_ds = mod_stats_grp.create_dataset('i_BMX', shape=(N_obj, 2), dtype=np.int, compression='gzip',
                                              compression_opts=4)  # Store which template is the bestmatch. [i_z, i_tpl]
        for i_obj in range(N_obj):
            bmx_flat = np.argmax(np.ravel(l[i_obj]))
            bmx_ds[i_obj, 0] = bmx_flat / N_tpl  # i_z bmx
            bmx_ds[i_obj, 1] = bmx_flat % N_tpl  # i_tpl bmx

        # 5.1.3 - Percentiles
        log.debug('Calculating Mass percentiles. t = %3.2f' % (time.time() - t0))
        p_grp.create_dataset('percentiles', shape=(N_obj,), data=np.array(
            [tuple(percentiles(np.ravel(s_Mass[i_obj]), np.ravel(likelihood_ds[i_obj]), perc)) for i_obj in range(N_obj)],
            dtype=perc_dt), compression='gzip', compression_opts=4)  # Percentiles dataset

        # 5.1.4 - Average
        log.debug('Calculating Mass average. t = %3.2f' % (time.time() - t0))
        p_grp.create_dataset('AVG', shape=(N_obj,), data=(s_Mass * likelihood_ds).sum(axis=2).sum(axis=1),
                             compression='gzip', compression_opts=4)

        # 5.2 - Library properties statistics
        for p_ in lib.properties.dtype.names:
            if lib.properties.dtype[p_].kind == 'f':  # Do statistics only for float properties.
                p_grp = lib_stats_grp.create_group(p_)  # Create a group for each prop

                # 5.2.1 - Check if property have NaN invalid elements.
                not_nan = np.isfinite(lib.properties[p_])
                are_nan = ~not_nan
                if are_nan.sum() > 0:
                    p_isnan = True
                    prop = np.ma.masked_array(lib.properties[p_], mask=are_nan)
                else:
                    p_isnan = False
                    prop = np.ma.masked_array(lib.properties[p_], mask=np.zeros_like(lib.properties[p_], dtype=np.bool))

                # 5.2.2 - Percentiles
                log.debug('Calculating property %s percentiles. t = %3.2f' % (p_, time.time() - t0))
                p_grp.create_dataset('percentiles', data=np.array(
                    [tuple(percentiles(prop.compressed(), l_T[i][not_nan], perc)) for i in range(N_obj)], dtype=perc_dt),
                                     compression='gzip', compression_opts=4)  # Percentiles dataset

                # 5.2.3 - Average
                log.debug('Calculating property %s average. t = %3.2f' % (p_, time.time() - t0))
                p_grp.create_dataset('AVG',
                                     data=np.sum(l_T * prop, axis=1) / np.sum(l_T[:, not_nan], axis=1), compression='gzip',
                                     compression_opts=4)

                # 5.2.4 - -999 probability (or NaN)
                if p_isnan:  # Will do this only if library property have NaNs.
                    log.debug('Calculating property %s NaN probability. t = %3.2f' % (p_, time.time() - t0))
                    p_grp.create_dataset('pNaN', data=np.sum(l_T[:, are_nan], axis=1), compression='gzip',
                                         compression_opts=4)

    f.close()

    log.info('Finished. Took %3.2f seconds.' % (time.time() - t0))
