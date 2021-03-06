'''
Created on Sep 18, 2012

@author: william

    This script converts STARLIGHT-SDSS input + output + table files
    to a hdf5 template library.  
'''

import re
import os
import time
import argparse
import logging
import h5py

import numpy as np

import atpy
from numpy.ma import mrecords
import pystarlight.io.starlighttable
import cosmocalc

from pystarlight.util.constants import Mpc_cm, L_sun

import magal.core.log

from magal.io.readfilterset import FilterSet
from magal.photometry.syntphot import photoconv
from magal.photometry.syntphot import spec2filterset
from magal.util.cosmo import zcor
from magal.util.matchs import matchobjs
from magal.io.hdf5util import inithdf5, read_filterhdf5

from magal.core.version import _magal_name_, _magal_version_


def argparser():
    
    ''' Defines the input and help options of the program... '''
    
    parser = argparse.ArgumentParser(description='Convert STARLIGHT-SDSS input + output + table output to a hdf5 library database.')
    
    parser.add_argument('-i', metavar='inputfilelist.txt', type=str, nargs=1,
                        help='Input file list with 2 columns: SDSS and Starlight files', required=True)    
    parser.add_argument('-so', metavar='/path/to/starlight_directory', type=str, nargs=1,
                        help='The directory where Starlight Output files are stored', required=True)
    parser.add_argument('-si', metavar='/path/to/sdss_directory', type=str, nargs=1,
                        help='The directory where SDSS Input .txt files are stored', required=True)
    parser.add_argument('-st', metavar='/path/to/tables_directory', type=str, nargs=1,
                        help='The directory where output Tables are stored', required=True)
    parser.add_argument('-b', metavar='BS', type=str, nargs=1, help='Base identifier (e.g. BS)', required=True)
    parser.add_argument('-d', metavar='DR?', type=str, nargs=1, help='SDSS Data Release (e.g. DR7)', required=True)
    parser.add_argument('-sa', metavar='sample', type=str, nargs=1, help='SDSS Sample (e.g. 926246)', required=True)
    parser.add_argument('-o', metavar='database.hdf5', type=str, nargs=1, help='Database output file', required=True)
    parser.add_argument('-f', metavar='curves.hdf5', type=str, nargs=1, help='Filter transmission curves file', required=True)
    parser.add_argument('-z_ini', metavar='0.0', type=float, nargs=1, help='Initial redshift', required=True)
    parser.add_argument('-z_fin', metavar='2.0', type=float, nargs=1, help='Final redshift', required=True)
    parser.add_argument('-dz', metavar='0.01', type=float, nargs=1, help='Delta redshift', required=True)
    parser.add_argument('-v', '-v', action='count')
    parser.add_argument('--version', action='version', version='%s version %s' % (_magal_name_, _magal_version_))
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    t0 = time.time() # To measure execution times.
    
    args = argparser()
    
    #logging
    from magal.core.log import setConsoleLevel
    setConsoleLevel(logging.DEBUG)
    
    log = logging.getLogger('magal.starlight2photdb')
    
    # 0 - Main definitions
    in_dir = args.si[0]
    out_dir = args.so[0]
    tables_dir = args.st[0]
    aux = np.loadtxt(args.i[0], dtype = np.str).T
    infiles = aux[0]
    outfiles = aux[1]
    el_file = '%s/sample.F%s.%s.f.lines.dat.BS.bz2' % (tables_dir, args.d[0], args.sa[0])
    syn01_file = '%s/sample.F%s.%s.f.Starlight.SYN01.tab.BS.bz2' % (tables_dir, args.d[0], args.sa[0])
    syn02_file = '%s/sample.F%s.%s.f.Starlight.SYN02.tab.BS.bz2' % (tables_dir, args.d[0], args.sa[0])
    syn03_file = '%s/sample.F%s.%s.f.Starlight.SYN03.tab.BS.bz2' % (tables_dir, args.d[0], args.sa[0])
    syn04_file = '%s/sample.F%s.%s.f.Starlight.SYN04.tab.BS.bz2' % (tables_dir, args.d[0], args.sa[0])
    db_file = args.o[0]
    filter_file = args.f[0]
    filterid = os.path.basename(filter_file).split('.')[0]
    z_from = np.float(args.z_ini[0])
    z_to = np.float(args.z_fin[0])
    z_step = np.float(args.dz[0])
    
    # 1 - Read files
    # 1.1 - Filter
    db_f = read_filterhdf5(filter_file)
    
    # 1.2 - Emission lines
    log.debug('Reading elines...')
    t_start = time.time()
    tb = atpy.Table(el_file, type='starlight_el')
    log.debug('Took %3.2f seconds.' % (time.time() - t_start))
    
    # 1.3 - SYN0[1-4]
    log.debug('Reading SYN files...')
    t_start = time.time()
    
    tsyn01 = atpy.Table(syn01_file, type='starlight_syn01', include_names = ('id', 'A_V', 'v0', 'vd', 'SN_w', 'SN_n'))
    tsyn02 = atpy.Table(syn02_file, type='starlight_syn02', include_names = ('id', 'at_flux', 'at_mass', 'aZ_flux', 'aZ_mass', 'am_flux', 'am_mass'))
    tsyn03 = atpy.Table(syn03_file, type='starlight_syn03', include_names = ('id', 'M2L_r'))
    tsyn04 = atpy.Table(syn04_file, type='starlight_syn04', include_names = ('id', 'Mcor_fib', 'Mcor_gal', 'DL_Mpc', 'Mini_fib', 'z'))

    log.debug('Took %3.2f seconds.' % (time.time() - t_start))
    
    # 1.4 - Join tables

    log.debug('Joining tables...')
    
    id_list = matchobjs(outfiles, tsyn04['id'])
    
    for key in tsyn01.keys():
        if(key != 'id'): tb.add_column(key, tsyn01[key])
    for key in tsyn02.keys():
        if(key != 'id'): tb.add_column(key, tsyn02[key])
    for key in tsyn03.keys():
        if(key != 'id'): tb.add_column(key, tsyn03[key])
    for key in tsyn04.keys():
        if(key != 'id' and key != 'DL_Mpc' and key != 'Mini_fib'): tb.add_column(key, tsyn04[key])
    
    # 0.1 - Init HDF5 file
    log.debug('Creating and populating DB file...')
    
    db = inithdf5(db_file)
        # Tables group
    db.create_group('/tables/')    
        # Filtersystem groups
    for filterid in db_f.keys():
        db.create_group('/%s/' % filterid)
            # CCD groups
        for ccd in db_f.get(filterid).keys():
            db.create_group('/%s/%s/' % (filterid, ccd))
        
    # 2 - Write tables to hdf5 file
    log.debug('\tTables...')
    properties = tb[id_list].copy()
    data = np.ma.array(properties.copy(),
                       mask=[tuple(np.append(False, np.less(list(p_)[1:], -998))) for p_ in properties],
                       fill_value=np.nan)
    # r = data.view(mrecords.mrecarray)
    ds_tb = db.create_dataset(name = '/tables/properties', data=data.filled())

    db.create_dataset(name = '/tables/z', data = np.arange(z_from, z_to, z_step) )

    
    # 3 - Photometry
    # 3.1 -
    log.debug('\tPhotometry...')
    for filterid in db_f.keys():
        for ccd in db_f.get(filterid).keys():
            # Read filtercurves
            f = FilterSet() # Init filterset object
            f.read(filter_file, path='/%s/%s' % (filterid, ccd))
            f.calc_filteravgwls()
            
            # db shape
            Nz = len(np.arange(z_from, z_to, z_step))
            Ngal = len(infiles)
            Nl = len(f.filteravgwls)
            
            # Library
            db.create_dataset(name = '/%s/%s/library' % (filterid, ccd), shape = (Nz,Ngal,Nl), dtype = np.dtype([('m_ab', np.float), ('e_ab', np.float)]) )
            db.create_dataset(name = '/%s/%s/filtercurves' % (filterid, ccd),  data = np.array(f.filterset, dtype=([('ID_filter', '|S32'), ('wl', '<f4'), ('transm', '<f4')])))
            aux = np.zeros(shape = len(f.filteravgwls), dtype = ([('ID_filter', '|S32'), ('wl_central', np.float)]))
            aux['ID_filter'] = np.unique(f.filterset['ID_filter'])
            aux['wl_central'] = f.filteravgwls
            db.create_dataset(name = '/%s/%s/filterset' % (filterid, ccd),  data = aux)
    
    # To convert units to L\odot / M\odot \AA:
    aux_units = (4 * np.pi * np.power(tsyn04['DL_Mpc'][id_list] * Mpc_cm,2)) / (np.power(10,tsyn04['Mini_fib'][id_list]) * L_sun)

    for i_file in range(Ngal):
        # Read starlight output
        model_file = '%s/%s' % (out_dir, outfiles[i_file])
        tm = atpy.TableSet(model_file, type='starlightv4')
        model_spec = np.copy(tm.spectra.data.view(dtype = np.dtype([('wl', tm.spectra.l_obs.dtype), ('f_obs', '<f8'), ('flux', tm.spectra.f_syn.dtype), ('f_wei', '<f8'), ('Best_f_SSP', '<f8')])))
        #     and input..
        obs_file = '%s/%s' % (in_dir, infiles[i_file])
        ts = atpy.TableSet(obs_file, type='starlight_input')
        obs_spec = np.copy(ts.starlight_input.data.view(dtype = np.dtype([('wl', '<f8'), ('flux', '<f8'), ('error', '<f8'), ('flag', '<i8')])))

        model_spec['flux'] = model_spec['flux'] * tm.keywords['fobs_norm'] * 1e-17 * aux_units[i_file]
        obs_spec['flux'] = obs_spec['flux'] * 1e-17 * aux_units[i_file]
        obs_spec['error'] = obs_spec['error'] * 1e-17 * aux_units[i_file]

        for ccd in db_f.get(filterid).keys():
            # Get the filterid filter.
            f.read(filter_file, path='/%s/%s' % (filterid, ccd))

            # For each defined redshift, eval the photometry and store on the database.
            db_m = db.get('/%s/%s/%s' % (filterid, ccd, 'library'))
            i_z = 0
            for z in np.arange(z_from, z_to, z_step):

                if z == 0:
                    d_L = 3.08567758e19 # 10 parsec in cm
                else:
                    d_L = cosmocalc.cosmocalc(z)['DL_cm']

                k_cosmo = L_sun  / ( 4 * np.pi * np.power(d_L,2) )

                O = zcor(obs_spec, z)
                O['flux'] = O['flux'] * k_cosmo
                O['error'] = O['error'] * k_cosmo

                M = zcor(model_spec, z)
                M['flux'] = M['flux'] * k_cosmo

                print 'i_file', i_file
                x = spec2filterset(f.filterset, O, M)

                db_m[i_z, i_file] = x
                i_z = i_z + 1
            

    db.close()
    
    log.debug('Took %3.2f seconds.' % (time.time() - t_start))