'''
Created on Dec 10, 2012

@author: william

    This script converts Base files to a PARAMETRIC hdf5 template library based on the given SFHs.  
'''

import time
import argparse
import matplotlib.pyplot as plt

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
setConsoleLevel(logging.CRITICAL)

def argparser():
    
    ''' Defines the input and help options of the program... '''
    
    parser = argparse.ArgumentParser(description='Convert STARLIGHT-SDSS input + output + table output to a hdf5 library database.')
    
    parser.add_argument('-B', metavar='BaseFile.txt', type=str, nargs=1,
                        help='Input base file given on STARLIGHT standards. See STARLIGHT manual for details.', required=True)    
    parser.add_argument('-Bd', metavar='/path/to/base_directory', type=str, nargs=1,
                        help='The directory where Base files are stored', required=True)
#    parser.add_argument('-si', metavar='/path/to/sdss_directory', type=str, nargs=1,
#                        help='The directory where SDSS Input .txt files are stored', required=True)
#    parser.add_argument('-st', metavar='/path/to/tables_directory', type=str, nargs=1,
#                        help='The directory where output Tables are stored', required=True)
#    parser.add_argument('-b', metavar='BS', type=str, nargs=1, help='Base identifier (e.g. BS)', required=True)
#    parser.add_argument('-d', metavar='DR?', type=str, nargs=1, help='SDSS Data Release (e.g. DR7)', required=True)
#    parser.add_argument('-sa', metavar='sample', type=str, nargs=1, help='SDSS Sample (e.g. 926246)', required=True)
    parser.add_argument('-o', metavar='database.hdf5', type=str, nargs=1, help='Database output file', required=True)
    parser.add_argument('-f', metavar='curves.hdf5', type=str, nargs=1, help='Filter transmission curves file', required=True)
    parser.add_argument('-z_ini', metavar='0.0', type=float, nargs=1, help='Initial redshift', required=True)
    parser.add_argument('-z_fin', metavar='2.0', type=float, nargs=1, help='Final redshift', required=True)
    parser.add_argument('-dz', metavar='0.01', type=float, nargs=1, help='Delta redshift', required=True)
    parser.add_argument('-v', '-v', action='count')
    parser.add_argument('--version', action='version', version='%s version %s' % (_magal_name_, _magal_version_))
    
    args = parser.parse_args()
    return args

#if __name__ == '__main__':

t0 = time.time() # To measure execution times.

#args = argparser()

#### Read Base File

#basedir = '/home/william/BasesDir/'
basedir = '/Users/william/mestrado/BasesDir/'
basefile = 'Base.bc03.Padova1994.chab.All'


bt = atpy.Table(basedir+basefile, basedir, read_basedir=True, type='starlightv4_base')

Z = np.unique(bt['Z_base']) 
ages = np.unique(bt['age_base'])


#### Read Filter Filter
filter_file = '/Users/william/doutorado/photo_filters/Alhambra_24.hdf5'
#filter_file = '/home/william/doutorado/photo_filters/sdss_ugriz.hdf5'
#filter_file = '/home/william/doutorado/photo_filters/Alhambra_24.hdf5'
db_f = read_filterhdf5(filter_file)

#################################################################
#
#        Library A = SSPs + A_V
#
#################################################################
#
#print 'Starting Library A...'
#
#z_from = 0.01
#z_to = 0.7
#z_step = 0.01
#
#av_from = 0.0
#av_to = 2.0
#av_step = .1
#av_sampling = np.arange(av_from, av_to, av_step)
#
##db_file = 'lib_csp_A.hdf5'
#db_file = '/home/william/lib_csp_B_sdss_z0p7.hdf5'
#try:
#    os.unlink(db_file)
#except:
#    pass
#db = inithdf5(db_file)
## Tables group
#db.create_group('/tables/')    
## Filtersystem groups
#for filterid in db_f.keys():
#    db.create_group('/%s/' % filterid)
## CCD groups
#for ccd in db_f.get(filterid).keys():
#    db.create_group('/%s/%s/' % (filterid, ccd))
#
#
## 1 - Tables
### 1.1 - redshift
#db.create_dataset(name = '/tables/z', data = np.arange(z_from, z_to, z_step) )
#
### 1.2 - properties
#aux_basesize = bt.shape[0]
#aux0 = drop_fields(bt.data.copy(), ['YA_V', 'aFe', 'l_ssp', 'f_ssp'])
#for i_av in range(len(av_sampling)):
#    aux1 = np.array(np.ones(len(bt['sspfile'])) * av_sampling[i_av], dtype=np.dtype([('AV', np.float)]))
#    aux2 = merge_arrays([aux0, aux1], flatten=True)
#    if i_av == av_sampling[0]:
#        aux_prop = db.create_dataset(name = '/tables/properties', shape=(aux_basesize * len(av_sampling),), dtype = aux2.dtype)
#    aux_prop[i_av * aux_basesize: (i_av+1) * aux_basesize] = aux2
#    
## 2 - Photometry
#for filterid in db_f.keys():
#    for ccd in db_f.get(filterid).keys():
#        
#        # 2.1 - Read filtercurves
#        f = FilterSet() # Init filterset object
#        f.read(filter_file, path='/%s/%s' % (filterid, ccd))
#        f.calc_filteravgwls()
#        
#        # db shape
#        Nz = len(np.arange(z_from, z_to, z_step))
#        Ngal = aux_basesize * len(av_sampling)
#        Nl = len(f.filteravgwls)
#        
#        # 2.2 - Create library
#        db.create_dataset(name = '/%s/%s/library' % (filterid, ccd), shape = (Nz,Ngal,Nl), dtype = np.dtype([('m_ab', np.float), ('e_ab', np.float)]) )
#        db.create_dataset(name = '/%s/%s/filtercurves' % (filterid, ccd),  data = np.array(f.filterset, dtype=([('ID_filter', '|S32'), ('wl', '<f4'), ('transm', '<f4')])))
#        aux = np.zeros(shape = len(f.filteravgwls), dtype = ([('ID_filter', '|S32'), ('wl_central', np.float)]))
#        aux['ID_filter'] = np.unique(f.filterset['ID_filter'])
#        aux['wl_central'] = f.filteravgwls
#        db.create_dataset(name = '/%s/%s/filterset' % (filterid, ccd),  data = aux)
#
## 2.3 - Populate library
#### Reddening Law == CCM ###
#q = calc_redlaw(bt[0]['l_ssp'], 3.1, 'CCM')
#
#
##aux_Base = drop_fields(bt.data, np.array(bt.data.dtype.names)[np.bitwise_and(np.array(bt.data.dtype.names) != 'l_ssp', np.array(bt.data.dtype.names) != 'f_ssp')])
#
#for i_base in range(aux_basesize):
#    print 'i_base', i_base
#    aux_Nl = len(bt.data[i_base]['l_ssp'])
#    aux_Bspec = np.zeros(aux_Nl, dtype = np.dtype([('wl', np.float), ('flux', np.float)]))
#    aux_Bspec['wl'] = bt[i_base]['l_ssp']
#    aux_Bspec['flux'] = bt[i_base]['f_ssp']
#    
#    i_av = 0
#    for av in av_sampling:
#        obs_spec = aux_Bspec.copy()
#        obs_spec['flux'] = (obs_spec['flux'] * 10**(-.4 * av * q))
#        
##        kk = np.bitwise_and(obs_spec['wl'] > 3000, obs_spec['wl'] < 9000)
#        i_z = 0
#        for ccd in db_f.get(filterid).keys():
#            # Get the filterid filter.
#            f.read(filter_file, path='/%s/%s' % (filterid, ccd))
#            
#            # For each defined redshift, eval the photometry and store on the database.
#            db_m = db.get('/%s/%s/%s' % (filterid, ccd, 'library'))
#            i_z = 0
#            for z in np.arange(z_from, z_to, z_step):
#                
#                if z == 0:
#                    d_L = 3.08567758e19 # 10 parsec in cm
#                else:
#                    d_L = cosmocalc.cosmocalc(z)['DL_cm']
#                    
#                k_cosmo = L_sun  / ( 4 * np.pi * np.power(d_L,2) ) 
#    
#                O = zcor(obs_spec, z)
#                O['flux'] = O['flux'] * k_cosmo
#
#                x = spec2filterset(f.filterset, O)
#                
##                print 'z', z
##                print aux0[i_base]
##                
##                plt.figure(1)
##                plt.clf()
##                plt.plot(O['wl'][kk], O['flux'][kk])
##                plt.figure(2)
##                plt.clf()
##                plt.plot(f.filteravgwls, x, '.')
#     
#                
#                db_m[i_z, i_base + aux_basesize*i_av] = x
#                i_z = i_z + 1
#        
#        i_av = i_av + 1
#                
##    print i_base
#        
#db.close()




#################################################################
#
#        Library B = 
#
#################################################################


# 0 - Library Parameters
libfile = '/Users/william/Downloads/databases/lib_csp_B_test.hdf5'

z_from = 0.005
z_to = .01 #
z_step = .01
Nz = len(np.arange(z_from, z_to, z_step))

libfile = '/Users/william/Downloads/databases/lib_csp_B_test_%3.2f_%3.2f.hdf5' % (z_step, z_to)

#From Bundy 2006 PhD thesis.
# The templates in this grid represent the assumed priors in this
# Bayesian technique and span 4 dimensions in parameter space: star
# formation history (\tau), age, metallicity, and dust content.

tau = np.array([0.04, 0.14, 0.40, 0.54, 0.60, 0.70, 0.71, 0.87, 0.92, 0.93, 1.39, 1.62,
                1.75, 1.75, 1.87, 1.99, 2.14, 2.63, 2.74, 3.06, 3.19, 4.08, 4.12, 4.79,
                5.26, 5.50, 6.09, 6.28, 6.98, 7.71, 7.88, 7.96, 8.17, 8.89, 9.10])[0:1] * 1e9  # For Gyr --> yr conversion

t0 = np.array([0.67, 0.98, 1.28, 2.74, 4.36, 5.06, 5.16, 6.15, 6.46, 6.75, 6.98, 7.97,
               8.71, 9.10, 9.53, 9.70])[0:1]  * 1e9  # For Gyr --> yr conversion

#Z = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]

N_av = 1
extinction = np.linspace(0,2, N_av)

lib_size = len(tau) * len(t0) * len(Z) * len(extinction) # \tau, age, metallicity and dust content.

print 'size> ', len(tau), len(t0), len(Z), len(extinction), lib_size

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
model_vec = {}
for ccd in ccds:
    db.create_group('/%s/%s/' % (filterid, ccd))
    f = FilterSet() # Init filterset object
    f.read(filter_file, path='/%s/%s' % (filterid, ccd))
    f.calc_filteravgwls()
    ccd_vec.update({ccd: f.filterset.copy()})
    Nl = len(f.filteravgwls)
    db.create_dataset(name = '/%s/%s/filtercurves' % (filterid, ccd),  data = np.array(f.filterset, dtype=([('ID_filter', '|S32'), ('wl', '<f4'), ('transm', '<f4')])))
    aux = np.zeros(shape = len(f.filteravgwls), dtype = ([('ID_filter', '|S32'), ('wl_central', np.float)]))
    aux['ID_filter'] = np.unique(f.filterset['ID_filter'])
    aux['wl_central'] = f.filteravgwls
    db.create_dataset(name = '/%s/%s/filterset' % (filterid, ccd),  data = aux)
    db_m = db.create_dataset('/%s/%s/%s' % (filterid, ccd, 'library'), shape = (Nz,lib_size, Nl), dtype = np.dtype([('m_ab', np.float), ('e_ab', np.float)]) )
    model_vec.update({ccd: np.empty(shape = (Nz,lib_size, Nl), dtype = np.dtype([('m_ab', np.float), ('e_ab', np.float)])) })




# 1 - Tables
## 1.1 - redshift
db.create_dataset(name = '/tables/z', data = np.arange(z_from, z_to, z_step) )

## 1.2 - properties
aux_dtype = np.dtype([('tau', np.float), ('t0', np.float), ('Z', np.float), ('A_V', np.float)])
prop = db.create_dataset(name = '/tables/properties', shape=(lib_size,), dtype = aux_dtype)
n_model = 0
for Z_lib in Z:
    for t_start in t0:
        for tau_eb in tau:
            for av in extinction:
                prop[n_model] = (tau_eb, t_start, Z_lib, av)
                n_model = n_model + 1


q = calc_redlaw(bt[Z == Z[0]]['l_ssp'][0], 3.1, 'CCM')
i_model = -1
aux_t = time.time()
aux_spec = np.zeros_like(bt[0]['f_ssp'], dtype = np.dtype([('wl', np.float), ('flux', np.float)]))
aux_spec['wl'] = bt[0]['l_ssp'][0]
bt_data = bt.data.copy()
for Z_lib in Z:
    print 'Z> ', Z_lib
    for t_start in t0:
        print 't_start> ', t_start
        for tau_eb in tau:
            print 'tau_eb> ', tau_eb
            for av in extinction:
                print 'av> ', av
                i_model += 1
                print 'i_model>', i_model
                print 't_model>', time.time() - aux_t
                aux_t = time.time()
                model = n_component(ages)
                model.add_exp(t_start, tau_eb, 1)
                sfh = model.get_sfh()
                spec = aux_spec.copy()
                dt = (model.ages_end - model.ages_start)
                for i in range(len(sfh)):
                    if sfh[i] > 0:
                        spec['flux'] = spec['flux'] + bt_data[bt_data['Z_base'] == Z_lib]['f_ssp'][i] * sfh[i] * 10**(-.4*av*q.copy()) * dt[i] 

                i_z = -1
                for z in np.arange(z_from, z_to, z_step):
                    i_z += 1
                    if z == 0:
                        d_L = 3.08567758e19 # 10 parsec in cm
                    else:
                        d_L = cosmocalc.cosmocalc(z)['DL_cm']
                        
                    k_cosmo = L_sun  / ( 4 * np.pi * np.power(d_L,2) )  #TODO: I can have this calculated beforehand
        
                    O = zcor(spec, z)
                    O['flux'] = O['flux'] * k_cosmo
    
                    if i_z % 20 == 0: print 'z>', z
    
                    for ccd in ccds:
                        fset = ccd_vec[ccd]
                        #db_m = db.get('/%s/%s/%s' % (filterid, ccd, 'library'))
                        x = spec2filterset(fset, O)
                        model_vec[ccd][i_z, i_model] = x

for ccd in ccds:
    db_m = db.get('/%s/%s/%s' % (filterid, ccd, 'library'))
    db_m = model_vec[ccd]
        
db.close()