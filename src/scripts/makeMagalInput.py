'''
Created on Jul 4, 2013

@author: william

This script reads catalog filenames and columns information from a config file and returns a magal .hdf5 file.  
'''

import sys
import ast, ConfigParser


import numpy as np
import h5py
from magal.io.readfilterset import FilterSet
from ConfigParser import NoOptionError


try:
    config_file = sys.argv[1]
except:
    print 'Usage: %s config_file' % sys.argv[0]
    sys.exit(1)
    
config = ConfigParser.ConfigParser()
config.read(sys.argv[1]) 

###
###  0 - Open filterfile to check if catalog is valid
###

Filter = FilterSet()
Filter.read(config.get('InputGeneral', 'filter_file'), path='%s/%s' % 
            (config.get('InputGeneral', 'filter_name'), ast.literal_eval(config.get('InputMagnitudes', 'ccds')+',')[0]))
Filter.calc_filteravgwls()

###
###  1 - DataType verification
###

# - Magnitudes
columns = ast.literal_eval(config.get('InputMagnitudes', 'columns'))
filter_ids = sorted(columns.keys()) # IMPORTANT! Always use a sorted version of the array!
input_dt = [(key, np.float) for key in filter_ids]
input_dt.extend([('d'+key, np.float) for key in filter_ids])
input_icols = [columns[key][0] for key in filter_ids]
input_icols.extend([columns[key][1] for key in filter_ids])

input_icols = input_icols - np.ones_like(input_icols)

# - Properties
properties_icols = [int(config.get('InputProperties', 'unique_id_column'))] # First column is always the unique identificator.
properties_dt = [('_id', np.int64)] # Unique id always haves _id name
# -- Float properties
columns = ast.literal_eval(config.get('InputProperties', 'general')) # IMPORTANT! Always use a sorted version of the array!
prop_ids = sorted(columns.keys())
properties_dt.extend([(key, np.float) for key in prop_ids])
properties_icols.extend([columns[key] for key in prop_ids])
# -- FLAG properties (They should be integers. Maybe I would have to check which are int and which are not.)
try:
    columns = ast.literal_eval(config.get('InputProperties', 'flags')) # IMPORTANT! Always use a sorted version of the array!
    prop_ids = sorted(columns.keys())
    properties_dt.extend([(key, np.int) for key in prop_ids])
    properties_icols.extend([columns[key] for key in prop_ids])
except NoOptionError:
    print 'The input file does not have any flag associated.'

properties_icols = properties_icols - np.ones_like(properties_icols) 

# - Check if the columns of the input file are the same as on filter file.
filters_filterfile = np.unique(Filter.filterset['ID_filter'])

missing_filters = [fid for fid in filter_ids if not fid in filters_filterfile]
if len(missing_filters) > 0:
    print 'Filter(s) %s does not exists on %s.' % (missing_filters, config.get('InputGeneral', 'filter_file'))
    sys.exit(0)
if len(filter_ids) < len(filters_filterfile):
    print 'There are less filters on %s than %s columns.' % (sys.argv[1], config.get('InputGeneral', 'filter_file'))
    sys.exit(0)
    
    #TODO: n_cat == n_ccds
    

Nfilters = len(filters_filterfile)

###
###  2 - Library init
###

# - Create .hdf5 input file.

db = h5py.File(config.get('InputGeneral', 'output_file'), 'w')

###
###  3 - Open ascii tables and read the data
###
catalog_fnames = ast.literal_eval(config.get('InputMagnitudes', 'cat_files')+',')

NgalCum = {}
i_cat = -1
for name in catalog_fnames:
    i_cat += 1
    aux_fname = '%s/%s' % (ast.literal_eval(config.get('InputMagnitudes', 'path_files')),
                           name)
    ccd = ast.literal_eval(config.get('InputMagnitudes', 'ccds')+',')[i_cat]
    mag_data = np.loadtxt(aux_fname, dtype=input_dt, usecols=input_icols)
    prop_data = np.loadtxt(aux_fname, dtype=properties_dt, usecols=properties_icols)

    Ngal = len(mag_data)
    
    ###
    ###  4 - Save redshift information (if exists)
    ###
    
    try:
        z_col = ast.literal_eval(config.get('InputProperties', 'redshift'))
    except:
        z_col = False
        
    if z_col:
        z_col = z_col.values()[0]
        redshift_data = np.loadtxt(aux_fname, dtype=np.float, usecols=(z_col,))
        db_redshift = db.get('/%s/%s/tables/z' % (config.get('InputGeneral', 'filter_name'), ccd))
        if db_redshift is None:
            db_redshift = db.create_dataset(name = '/%s/%s/tables/z' % ((config.get('InputGeneral', 'filter_name'), ccd)),
                                  data = redshift_data,  
                                  dtype = np.float,
                                  chunks = True,
                                  maxshape = (None,) )
        else:
            db_redshift[NgalCum[ccd]:NgalCum[ccd]+Ngal] = redshift_data
            
        
        

    ###
    ###  5 - Save properties table data.
    ###
        
    db_properties = db.get('/%s/%s/tables/properties' % (config.get('InputGeneral', 'filter_name'), ccd))
    if db_properties is None:
        db_properties = db.create_dataset(name = '/%s/%s/tables/properties' % ((config.get('InputGeneral', 'filter_name'), ccd)),
                                          data = prop_data,  
                                          dtype = properties_dt,
                                          chunks = True,
                                          maxshape = (None,) )
    else:
        db_properties.resize((NgalCum[ccd]+Ngal,))
        db_properties[NgalCum[ccd]:NgalCum[ccd]+Ngal] = prop_data
    
    ###
    ###  6 - Load magnitudes to the table.
    ###
    
    db_data = db.get('/%s/%s/data' % (config.get('InputGeneral', 'filter_name'), ccd))
    if db_data is None:
        db_data = db.create_dataset(name = '/%s/%s/data' % ((config.get('InputGeneral', 'filter_name'), ccd)),
                                    shape = (Ngal,Nfilters),
                                    dtype = np.dtype([('m_ab', np.float), ('e_ab', np.float)]),
                                    chunks = True,
                                    maxshape = (None, Nfilters) )
        NgalCum[ccd] = 0
    else:
        db_data.resize(((NgalCum[ccd]+Ngal,Nfilters)))
    
    for i_obj in range(Ngal):
        x = np.zeros(Nfilters, dtype = np.dtype([('m_ab', '<f8'), ('e_ab', '<f8')])) 
        x['m_ab'] = [mag_data[i_obj][fil] for fil in filters_filterfile]  
        x['e_ab'] = [mag_data[i_obj]['d%s' % fil] for fil in filters_filterfile] 
        badpxl = np.bitwise_or(x['m_ab'] >= 97, x['m_ab'] <= -97) #FIXME: Should be passed by the config file.
        x['m_ab'][badpxl] = np.inf
        x['e_ab'][badpxl] = np.inf
     
        db_data[i_obj+NgalCum[ccd]] = x
    
#         if i_obj % 100 == 0: print 'DEBUG: I\'m at i_obj, i_z --> %s' % i_obj

    NgalCum[ccd] += Ngal

db.close()