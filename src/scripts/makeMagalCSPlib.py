'''
Created on Dec 10, 2012

@author: william

    This script converts Base files to a PARAMETRIC hdf5 template library based on the given SFHs.  
'''

import time
import argparse
import ConfigParser
# import matplotlib.pyplot as plt

from magal.core.version import _magal_name_, _magal_version_

import atpy
import numpy as np
import pystarlight.io #@UnusedImport

from magal.util.stellarpop import n_component
from pystarlight.util.redenninglaws import calc_redlaw
from magal.io.hdf5util import inithdf5, read_filterhdf5
from numpy.lib.recfunctions import drop_fields, merge_arrays, rename_fields
from magal.io.readfilterset import FilterSet
import cosmocalc
from pystarlight.util.constants import L_sun
from magal.util.cosmo import zcor
from magal.photometry.syntphot import spec2filterset
import os

#logging
import logging
from magal.core.log import setConsoleLevel
import itertools
import sys
import ast
setConsoleLevel(logging.CRITICAL)


#if __name__ == '__main__':

t0 = time.time() # To measure execution times.

config = ConfigParser.ConfigParser()
config.read(sys.argv[1])
#TODO: Put some checking/error handling here.

#### Read Base File
basedir = config.get('LibraryGeneral', 'bases_dir')
basefile = config.get('LibraryGeneral', 'base_file')

bt = atpy.Table(basedir+basefile, basedir, read_basedir=True, type='starlightv4_base')

Z = np.unique(bt['Z_base']) 
ages = np.unique(bt['age_base'])
# Check if there is a restriction on the metallicity passed by command line argument.
if len(sys.argv) > 3:
    Z = [Z[int(sys.argv[2])]]
    libfile = '%s_%s' % (config.get('LibraryGeneral', 'lib_file'), sys.argv[2])
else:
    libfile = config.get('LibraryGeneral', 'lib_file')

#### Read Filter Filter
filter_file = config.get('LibraryGeneral', 'filter_file')
db_f = read_filterhdf5(filter_file)

# 0 - Library Parameters
z_from = np.float(config.get('LibraryParameters', 'z_from'))
z_to = np.float(config.get('LibraryParameters', 'z_to'))
z_step = np.float(config.get('LibraryParameters', 'z_step'))

z_range = np.arange(z_from, z_to, z_step)
Nz = len(z_range)



#From Bundy 2006 PhD thesis.
# The templates in this grid represent the assumed priors in this
# Bayesian technique and span 4 dimensions in parameter space: star
# formation history (\tau), age, metallicity, and dust content.

tau = ast.literal_eval(config.get('LibraryParameters', 'tau'))
t0 = ast.literal_eval(config.get('LibraryParameters', 't0'))

av_range = ast.literal_eval(config.get('LibraryParameters', 'a_v'))

lib_size = len(tau) * len(t0) * len(Z) * len(av_range) # \tau, age, metallicity and dust content.

try:
    os.unlink(libfile)
except:
    pass
db = inithdf5(libfile)
# Tables group
db.create_group('/tables/')    
# Filtersystem groups
for filterid in db_f.keys():
    db.create_group('/%s/' % filterid)
# CCD groups
ccds = db_f.get(filterid).keys()
ccd_vec = {}
db_vec = {}
for ccd in ccds:
    db.create_group('/%s/%s/' % (filterid, ccd))
    f = FilterSet() # Init filterset object
    f.read(filter_file, path='/%s/%s' % (filterid, ccd))
    f.calc_filteravgwls()
    ccd_vec.update({ccd: f.filterset})
    Nl = len(f.filteravgwls)
    db.create_dataset(name = '/%s/%s/filtercurves' % (filterid, ccd),  data = np.array(f.filterset, dtype=([('ID_filter', '|S32'), ('wl', '<f4'), ('transm', '<f4')])))
    aux = np.zeros(shape = len(f.filteravgwls), dtype = ([('ID_filter', '|S32'), ('wl_central', np.float)]))
    aux['ID_filter'] = np.unique(f.filterset['ID_filter'])
    aux['wl_central'] = f.filteravgwls
    db.create_dataset(name = '/%s/%s/filterset' % (filterid, ccd),  data = aux)
    chunk_z = Nz #, #np.int(64*1024/(8*Nl))
    db_m = db.create_dataset('/%s/%s/%s' % (filterid, ccd, 'library'), shape = (Nz,lib_size, Nl), chunks = (chunk_z,1,Nl), compression=4, dtype = np.dtype([('m_ab', '<f4'), ('e_ab', '<f4')]) )
    db_vec.update({ccd: db_m})




# 1 - Tables
## 1.1 - redshift
db.create_dataset(name = '/tables/z', data = z_range )

## 1.2 - properties
aux_dtype = np.dtype([('tau', np.float), ('t0', np.float), ('Z', np.float), ('A_V', np.float)])
prop = db.create_dataset(name = '/tables/properties', shape=(lib_size,), dtype = aux_dtype)
n_model = 0

for Z_lib in Z:
    for t_start in t0:
        for tau_eb in tau:
            for av in av_range:
                prop[n_model] = (tau_eb, t_start, Z_lib, av)
                n_model = n_model + 1

def spec2filter_z(args):
    spec, ccd_filter, z = args
    if z == 0:
        d_L = 3.08567758e19 # 10 parsec in cm
    else:
        d_L = cosmocalc.cosmocalc(z)['DL_cm']
        
    k_cosmo = L_sun  / ( 4 * np.pi * np.power(d_L,2) ) 

    O = zcor(spec, z)
    O['flux'] = O['flux'] * k_cosmo

    x = spec2filterset(ccd_filter, O)
    
    return x

from multiprocessing import Pool
q = calc_redlaw(bt[Z == Z[0]]['l_ssp'][0], 3.1, 'CCM')
i_model = -1
aux_t = time.time()

for Z_lib in Z:
    for t_start in t0:
        print 't_start> ', t_start
        for tau_eb in tau:
            print 'tau_eb> ', tau_eb
            for av in av_range:
                print 'av> ', av
                i_model += 1
                print 'i_model> %i of %i' % (i_model, lib_size)
                print 't_model>', time.time() - aux_t
                aux_t = time.time()
                model = n_component(ages)
                model.add_exp(t_start, tau_eb, 1)
                sfh = model.get_sfh()
                spec = np.zeros(shape=np.shape(bt[0]['f_ssp']), dtype = np.dtype([('wl', np.float), ('flux', np.float)]))
                spec['wl'] = bt[bt['Z_base'] == Z_lib]['l_ssp'][0]
                dt = (model.ages_end - model.ages_start)
                for i in range(len(sfh)):
                    if sfh[i] > 0: spec['flux'] = spec['flux'] + bt[bt['Z_base'] == Z_lib][i]['f_ssp'] * sfh[i] * 10**(-.4*av*q.copy()) * dt[i] / bt[i]['Mstars']  #ADDED THE MSTARS TERM!!! a ver que pasa chiquillo... :)
                for ccd in ccds:
                    args = []
                    for z in z_range:
                        args.append((spec, ccd_vec[ccd], z))
                    pool = Pool()
                    result = pool.map(spec2filter_z, args)
                    pool.close()
                    pool.join()
                    db_vec[ccd][:,i_model,:] = result
                        
        
db.close()