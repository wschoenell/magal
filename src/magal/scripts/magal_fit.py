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
import sys
import os
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import ConfigParser
from ConfigParser import NoSectionError, NoOptionError
import ast
import multiprocessing

import h5py
import numpy as np
import time

from magal.core.log import logger
from magal.core.version import _magal_version_, _magal_updated_
from magal.core.exceptions import MAGALCLIError
from magal.io.readlibrary import Library
from magal.fit.stats import percentiles
from magal.io.readinput import Input
from magal.util.matchs import get_zslice
from magal.fit.magal_aux import chi2_parameters, chi2_wrapper

DEBUG = 0
TESTRUN = 0
PROFILE = 0


def md5_for_file(f, block_size=2 ** 20):
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()


def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    log = logging.getLogger('magal.main')

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % _magal_version_
    program_build_date = str(_magal_updated_)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

USAGE
''' % (program_shortdesc)

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)

        parser.add_argument("-i", "--inputfile", dest="input", help="Input configuration", required=True)
        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                            help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()
        if args.verbose > 0:
            log.setLevel('DEBUG')

        #### ALL the configfile variables are loaded here ####
        config = ConfigParser.ConfigParser()
        config.read(args.input)

        mag_file = os.path.expandvars(config.get('FitGeneral', 'mag_file'))  # Input File
        lib_file = os.path.expandvars(config.get('FitGeneral', 'lib_file'))  # Template ibrary File
        output_file = os.path.expandvars(config.get('FitGeneral', 'output_file'))  # Output File
        filter_sys = config.get('FitGeneral', 'filter_sys')  # Filtersystem (e.g. SDSS)
        ccd = config.get('FitGeneral', 'ccd')  # CCD number

        #If allow_overwrite is defined, we can overwrite the outputfile.
        try:
            allow_overwrite = config.getboolean('FitGeneral', 'allow_overwrite')
        except NoOptionError:
            allow_overwrite = False
        #Check if outputfile exists and act as defined by config.
        if os.path.exists(output_file):
            if allow_overwrite:
                try:
                    os.unlink(output_file)
                except:
                    print u'Cannot DELETE file {0:s}.'.format(output_file)
                    return 1
            else:
                print u'File {0:s} already exists.'.format(output_file)
                return 1

        # Look if there is a filter list to include.
        try:
            filters_include = ast.literal_eval(config.get('FitGeneral', 'filters_include'))
        except NoOptionError:
            log.debug('Fitting for all the filters.')
            filters_include = None

        # Look if there is a maximum number of objects to run. This is useful for testing.
        try:
            Nobj_max = config.getint('FitGeneral', 'Nobj_max')
        except NoOptionError:
            Nobj_max = 0
            log.info('Running magal fit to ALL objects.')

        # Check if the user request a simulation. For a simulation, the inputfile will be a Library template file.
        try:
            is_simulation = config.getboolean('FitSimulation', 'is_simulation')  # Is this a simulation?
            if is_simulation:
                try:
                    obj_z = config.getfloat('FitSimulation', 'obj_z')
                    print('Using z = %s from inputfile on this simulation.' % config.getfloat('FitSimulation', 'obj_z'))
                except NoOptionError:
                    print('Using z from inputfile tables.')
                try:
                    mass_field = config.get('FitSimulation', 'mass_field')
                except NoOptionError:
                    print('When running a simulation, you must specify a table field for the mass.')
                    sys.exit(1)
            else:
                mass_field = None
        except NoSectionError:
            is_simulation = False  # If there is no [FitSimulation] on config, just ignore this.

        try:
            Nz = config.getboolean('FitGeneral', 'Nz')  # If True, will fit at z = const.
            # Redshift will be get from inputfile property z.
        except NoOptionError:
            if is_simulation:
                Nz = False
            else:
                log.info('Evaluating for all redshift space!')
        ########

    except KeyboardInterrupt:
        print 'CTRL+C pressed... exiting...'  #TODO: move this to the logger.
        return 0

    except Exception, e:
        if DEBUG or TESTRUN:
            raise (e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

    # 1 - Load files
    # 1.1 - Inputfile
    if is_simulation:  # If this is a simulation, a Library will be the inputfile.
        inp = Library(mag_file)
        inp.get_filtersys(filter_sys, ccd)
        magal_type = 'simulation'
    else:
        inp = Input(mag_file)
        try:
            inp.get_filtersys(filter_sys, ccd)
        except KeyError:
            raise MAGALCLIError('Data not found! Are you sure that you do not want to Simulate?')
        magal_type = 'data'

    # 1.2 - Library
    lib = Library(lib_file)
    lib.get_filtersys(filter_sys, ccd)

    # 1.3 - Get objectlist
    if is_simulation:
        o_list = get_zslice(inp, obj_z)  # Need to get selected redshift from the library
    else:
        o_list = inp.data

    if Nobj_max > 0:
        o_list = o_list[:Nobj_max]
        if args.verbose:
            print 'Running magal fit to ONLY %i objects.' % Nobj_max
    elif args.verbose:
        print 'Running magal fit to ALL objects.'

    N_obj = len(o_list)

    if is_simulation:
        if mass_field is not None:
            mass_field = inp.properties[mass_field][:N_obj]  # Multiply by the mass_field.

    # 2 - Output file
    # 2.1 - init 
    try:
        f = h5py.File(output_file, mode='w-')
    except IOError:
        raise MAGALCLIError(u'File {0:s} already exists or could not be created.'.format(output_file))

    # 2.2 - Define some auxiliar data...
    f.attrs.create('ifile', mag_file)
    fp = open(mag_file, 'rb')
    f.attrs.create('ifile_md5', md5_for_file(fp))  # Store Inputfile MD5SUM
    f.attrs.create('lib', lib_file)
    fp = open(lib_file, 'rb')
    f.attrs.create('lib_md5', md5_for_file(fp))  # Store Library MD5SUM
    f.attrs.create('type', magal_type)
    f.attrs.create('version', '%s - %s' % (program_name, program_version))

    # Write full .ini file on model outputfile.
    aux_ini = open(args.input, 'r').read()
    ds = f.create_dataset('/ini_file', shape=(1,), dtype=h5py.new_vlen(str))
    ds[:] = aux_ini
    ### 

    # 2.3 - Data matrices
    # 2.3.1 - Shape is defined by (N_obj, N_z, N_lib)
    if Nz:
        aux_shape = (N_obj, 1, lib.library.shape[1])
    else:
        aux_shape = (N_obj, lib.library.shape[0], lib.library.shape[1])

    # 2.3.2 - Create datasets
    n_ds = f.create_dataset('%s/n' % (lib.path), shape=aux_shape, dtype=np.int, compression='gzip',
                            compression_opts=4)  # Number of used pixels on \chi2 calc
    s_ds = f.create_dataset('%s/s' % (lib.path), shape=aux_shape, compression='gzip',
                            compression_opts=4)  # Scaling-factor
    chi2_ds = f.create_dataset('%s/chi2' % (lib.path), aux_shape, compression='gzip', compression_opts=4)  # \Chi2
    mod_stats_grp = f.create_group(
        '%s/statistics/model' % (lib.path))  # Model likelihood statistics (i.e. average, percentiles, etc..)
    lib_stats_grp = f.create_group('%s/statistics/library' % (lib.path))  # Library likelihood statistics

    # 3 - Fit
    if Nz:
        N_z, N_tpl = (1, lib.library.shape[1])  # If the redshift comes from inputfile, N_z = 1.
        inp_z = inp.z
    else:
        N_z, N_tpl = lib.library.shape[:2]
        inp_z = lib.z

    if filters_include is not None:
        filterset_mask = np.sum([k == lib.filterset['ID_filter'] for k in filters_include], axis=0, dtype=np.bool)
    else:
        filterset_mask = np.ones(len(lib.filterset), dtype=np.bool)

    # 3.3 - Eval \chi^2
    params = chi2_parameters(Nz, is_simulation, o_list, N_obj, N_z, N_tpl, inp_z, filterset_mask, lib_file, filter_sys,
                             ccd, mass_field)
    pool = multiprocessing.Pool()
    print 'started pool'
    for result in pool.map(chi2_wrapper, params):
        i_obj, aux_n, aux_s, aux_chi2 = result
        n_ds[i_obj], s_ds[i_obj], chi2_ds[i_obj] = aux_n, aux_s, aux_chi2

    # 4 - From \chi^2, eval likelihood
    # 4.1 - Likelihood
    log.debug('Calculating Likelihood')
    l = np.exp(-0.5 * ( np.subtract(chi2_ds, np.min(chi2_ds, axis=2).reshape(N_obj, N_z, 1)) ))

    # 4.2 - Normalization
    l /= np.sum(l, axis=2).reshape(N_obj, N_z, 1)  #FIXME: Double sum here!
    likelihood_ds = f.create_dataset('%s/likelihood' % (lib.path), data=l, compression='gzip', compression_opts=4)

    l_T = np.sum(l, axis=1)  # Template likelihood
    norm = np.sum(l_T, axis=1)  # Template likelihood normalization
    l_T = np.array([l_T[i_obj] / norm[i_obj] for i_obj in range(N_obj)])

    f.create_dataset('%s/likelihood_template' % lib.path, data=l_T, compression='gzip', compression_opts=4)
    if N_z > 1:  # If we fit over the redshift space, calulate the redshift-compressed likelihood
        l_z = np.sum(l, axis=2)
        norm = np.sum(l_z, axis=1)
        l_z = np.array([l_z[i_obj] / norm[i_obj] for i_obj in range(N_obj)])
        lz_ds = f.create_dataset('%s/likelihood_redshift' % (lib.path), data=np.array(l_z), compression='gzip',
                                 compression_opts=4)  # @UnusedVariable

    # 5 - From likelihood, eval statistics...
    # 5.0 - Some definitions...
    perc = [2.5, 16, 50, 84, 97.5]  # Percentiles which we calculate.
    perc_dt = np.dtype([('%s' % p, np.float) for p in perc])

    # 5.1 - Model statistics
    # 5.1.1 - Scaling-factor "s".
    log.debug('Calculating scaling-factor statistics')
    p_grp = mod_stats_grp.create_group('s_Mass')  # Scale prop group
    s_Mass = np.divide(s_ds.value.copy(), -2.5)

    # 5.1.2 - Best match
    log.debug('Calculating Best Matches')
    bmx_ds = mod_stats_grp.create_dataset('i_BMX', shape=(N_obj, 2), dtype=np.int, compression='gzip',
                                          compression_opts=4)  # Store which template is the bestmatch. [i_z, i_tpl]
    for i_obj in range(N_obj):
        bmx_flat = np.argmax(np.ravel(l[i_obj]))
        bmx_ds[i_obj, 0] = bmx_flat / N_tpl  # i_z bmx
        bmx_ds[i_obj, 1] = bmx_flat % N_tpl  # i_tpl bmx

    # 5.1.3 - Percentiles
    log.debug('Calculating Mass percentiles')
    p_grp.create_dataset('percentiles', shape=(N_obj,), data=np.array(
        [tuple(percentiles(np.ravel(s_Mass[i_obj]), np.ravel(likelihood_ds[i_obj]), perc)) for i_obj in range(N_obj)],
        dtype=perc_dt), compression='gzip', compression_opts=4)  # Percentiles dataset

    # 5.1.4 - Average
    log.debug('Calculating Mass average')
    p_grp.create_dataset('AVG', shape=(N_obj,),
                         data=[np.sum(s_Mass[i_obj] * likelihood_ds[i_obj]) for i in range(N_obj)], compression='gzip',
                         compression_opts=4)

    # 5.2 - Library properties statistics
    for p_ in lib.properties.dtype.names:
        if lib.properties.dtype[p_].kind == 'f':  # Do statistics only for float properties.
            p_grp = lib_stats_grp.create_group(p_)  # Create a group for each prop

            # 5.2.1 - Check if property have NaN invalid elements.
            are_nan = np.isnan(lib.properties[p_])
            not_nan = np.invert(are_nan)
            if are_nan.sum() > 0:
                p_isnan = True
                prop = np.ma.masked_array(lib.properties[p_], mask=are_nan)
            else:
                p_isnan = False
                prop = np.ma.masked_array(lib.properties[p_], mask=np.zeros_like(lib.properties[p_], dtype=np.bool))

            # 5.2.2 - Percentiles
            log.debug('Calculating property %s percentiles' % p_)
            p_grp.create_dataset('percentiles', data=np.array(
                [tuple(percentiles(prop.compressed(), l_T[i][not_nan], perc)) for i in range(N_obj)], dtype=perc_dt),
                                 compression='gzip', compression_opts=4)  # Percentiles dataset

            # 5.2.3 - Average
            log.debug('Calculating property %s average' % p_)
            p_grp.create_dataset('AVG',
                                 data=np.sum(l_T * prop, axis=1) / np.sum(l_T[:, not_nan], axis=1), compression='gzip',
                                 compression_opts=4)

            # 5.2.4 - -999 probability (or NaN)
            if p_isnan:  # Will do this only if library property have NaNs.
                log.debug('Calculating property %s NaN probability' % prop)
                p_grp.create_dataset('pNaN', data=np.sum(l_T[:, are_nan], axis=1), compression='gzip',
                                     compression_opts=4)

    f.close()

    log.debug('Finished.')


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-v")
    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'scripts.magal.profile'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())