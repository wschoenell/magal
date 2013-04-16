'''
Created on Sep 18, 2012

@author: william

    This script converts STARLIGHT-SDSS input + output + table files
    to a hdf5 template library.  
'''

import os
import time
import argparse
import logging

import numpy as np

import atpy

import magal.core.log

from magal.io.readfilterset import FilterSet
from magal.photometry.syntphot import photoconv
from magal.io.hdf5util import inithdf5, read_filterhdf5

from magal.core.version import _magal_name_, _magal_version_


def argparser():
    
    ''' Defines the input and help options of the program... '''
    
    parser = argparse.ArgumentParser(description='Convert STARLIGHT-SDSS input + output + table output to a hdf5 library database.')
    
    parser.add_argument('-i', metavar='inputfilelist.txt', type=str, nargs=1,
                        help='Input file list with SDSS .fits files on each row', required=True)
    parser.add_argument('-o', metavar='database.hdf5', type=str, nargs=1, help='Database output file', required=True)
    parser.add_argument('-f', metavar='curves.hdf5', type=str, nargs=1, help='Filter transmission curves file', required=True)
#     parser.add_argument('-z_ini', metavar='0.0', type=float, nargs=1, help='Initial redshift', required=False)
#     parser.add_argument('-z_fin', metavar='2.0', type=float, nargs=1, help='Final redshift', required=False)
#     parser.add_argument('-dz', metavar='0.01', type=float, nargs=1, help='Delta redshift', required=False)
    parser.add_argument('-v', '-v', action='count')
    parser.add_argument('--version', action='version', version='%s version %s' % (_magal_name_, _magal_version_))
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    
    t0 = time.time() # To measure execution times.
    
    args = argparser()
    
    #logging
    from magal.core.log import setConsoleLevel
    
    log = logging.getLogger('magal.sdss2photdb')
    
    if args.v > 0: setConsoleLevel('DEBUG')
    else: setConsoleLevel('CRITICAL')
    
    # 0 - Main definitions
    t_start = time.time()
    infiles = np.loadtxt(args.i[0], dtype = np.str) #np.dtype((['filename', np.str]))
    db_file = args.o[0]
    filter_file = args.f[0]
    filterid = os.path.basename(filter_file).split('.')[0]
#    z_from = np.float(args.z_ini[0])
#    z_to = np.float(args.z_fin[0])
#    z_step = np.float(args.dz[0])
    
    # 1 - Read files
    # 1.1 - Filter
    db_f = read_filterhdf5(filter_file)
        
    # 1.1 - Init HDF5 file
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
    ds_prop = db.create_dataset(name = '/tables/properties', shape=(len(infiles), 1), dtype=np.dtype([('SPECOBJID', '|S22'), ('SOURCETYPE', '|S19'), ('Z', '>f4'), ('Z_ERR', '>f4')]))
#    db.create_dataset(name = '/tables/z', data = np.arange(z_from, z_to, z_step) )

    
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
            Nz = 1
            Nobj = len(infiles)
            Nl = len(f.filteravgwls)
            
            # Library
            db.create_dataset(name = '/%s/%s/library' % (filterid, ccd), shape = (Nz,Nobj,Nl), dtype = np.dtype([('m_ab', np.float), ('e_ab', np.float)]) )
            db.create_dataset(name = '/%s/%s/filtercurves' % (filterid, ccd),  data = np.array(f.filterset, dtype=([('ID_filter', '|S32'), ('wl', '<f4'), ('transm', '<f4')])))
            aux = np.zeros(shape = len(f.filteravgwls), dtype = ([('ID_filter', '|S32'), ('wl_central', np.float)]))
            aux['ID_filter'] = np.unique(f.filterset['ID_filter'])
            aux['wl_central'] = f.filteravgwls
            db.create_dataset(name = '/%s/%s/filterset' % (filterid, ccd),  data = aux)
    
    for i_file in range(Nobj):
            
        for ccd in db_f.get(filterid).keys():
            # Get the filterid filter.
            f.read(filter_file, path='/%s/%s' % (filterid, ccd))
            
            db_m = db.get('/%s/%s/%s' % (filterid, ccd, 'library'))
            i_z = 0
            
            x = photoconv()
            print infiles[i_file]
            fits = atpy.Table(infiles[i_file], type='fits', hdu='COADD')
            properties = atpy.Table(infiles[i_file], type='fits', hdu='SPECOBJ')
            ds_prop[i_file] = properties[['SPECOBJID', 'SOURCETYPE', 'Z', 'Z_ERR']]
                
            x = x.fromSDSSfits(f.filterset, fits)
            
            db_m[i_z, i_file] = x
            

    db.close()
    
    log.debug('Took %3.2f seconds.' % (time.time() - t_start))