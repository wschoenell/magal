#!/usr/bin/env python
'''
Created on Jul 4, 2013

@author: william

This script reads catalog filenames and columns information from a config file and returns a magal .hdf5 file.  
'''

import sys
import ast
import ConfigParser
from ConfigParser import NoOptionError

import numpy as np
import h5py

from magal.core.log import logger


log = logger(__name__)
log.setLevel('DEBUG')


def main():
    if len(sys.argv) < 2:
        print 'Usage: %s configuration.ini' % sys.argv[0]
        return 1

    # #### Config file #####
    try:
        config = ConfigParser.ConfigParser()
        config.read(sys.argv[1])
    except:
        print 'Usage: %s config_file' % sys.argv[0]
        sys.exit(1)

    # # General
    output_file = config.get('InputGeneral', 'output_file')
    filterset_name = config.get('InputGeneral', 'filterset_name')

    ## Magnitudes
    path_files = ast.literal_eval(config.get('InputMagnitudes', 'path_files'))
    catalog_fnames = ast.literal_eval(config.get('InputMagnitudes', 'cat_files'))
    ccds = ast.literal_eval(config.get('InputMagnitudes', 'ccds'))

    try:
        catalog_delimiter = ast.literal_eval(config.get('InputGeneral', 'delimiter'))
    except:
        catalog_delimiter = None

    columns = ast.literal_eval(config.get('InputMagnitudes', 'columns'))
    filter_ids = sorted(columns.keys())  # IMPORTANT! Always use a sorted version of the array!
    Nfilters = len(filter_ids)
    input_dt = [(key, np.float) for key in filter_ids]
    input_dt.extend([('d' + key, np.float) for key in filter_ids])
    input_icols = [columns[key][0] for key in filter_ids]
    input_icols.extend([columns[key][1] for key in filter_ids])

    input_icols = input_icols - np.ones_like(input_icols)


    ### Galactic extinction correction.
    # TODO: implement the option to use ra,dec or l,b and Schlegel maps of obstools.
    try:
        extinction = ast.literal_eval(config.get('InputGeneral', 'extinction'))
        extinction_dt = [(key, np.float) for key in extinction.keys()]
        extinction_icols = [extinction[key] for key in extinction.keys()]
        extinction_icols = extinction_icols - np.ones_like(extinction_icols)
        if extinction.keys() != columns.keys():
            raise 'Extinction and data columns are different!'
    except NoOptionError:
        extinction = False
        print 'The input file does not have Galactic extinction correction to the magnitudes.'

    ## Properties
    unique_id_column = config.getint('InputProperties', 'unique_id_column')
    properties_icols = [unique_id_column]  # First column is always the unique identificator.
    properties_dt = [('_id', np.int64)]  # Unique id always haves _id name

    ### Float properties
    columns = ast.literal_eval(config.get('InputProperties', 'general'))
    prop_ids = sorted(columns.keys())  # IMPORTANT! Always use a sorted version of the array!
    properties_dt.extend([(key, np.float) for key in prop_ids])
    properties_icols.extend([columns[key] for key in prop_ids])

    ### FLAG properties (They should be integers. Maybe I would have to check which are int and which are not.)
    try:
        columns = ast.literal_eval(config.get('InputProperties', 'flags'))
        prop_ids = sorted(columns.keys())  # IMPORTANT! Always use a sorted version of the array!
        properties_dt.extend([(key, np.int) for key in prop_ids])
        properties_icols.extend([columns[key] for key in prop_ids])
    except NoOptionError:
        print 'The input file does not have any flag associated.'
    properties_icols = properties_icols - np.ones_like(properties_icols)

    ### Redshift
    try:
        z_col = config.getint('InputProperties', 'redshift')
        z_col -= 1
    except NoOptionError:
        z_col = False

    ##### --- #####
    ###
    ###  0 - Library init
    ###

    # - Create .hdf5 input file.
    db = h5py.File(output_file, 'w')
    # - Write full .ini file on inputfile.
    aux_ini = open(sys.argv[1], 'r').read()
    str_type = h5py.new_vlen(str)
    ds = db.create_dataset('/ini_file', shape=(1,), dtype=str_type)
    ds[:] = aux_ini
    ## - Define few file attributes
    db.attrs.create('filters', filter_ids)

    ###
    ###  1 - Open ascii tables and read the data
    ###
    NgalCum = {}
    i_cat = -1
    for name in catalog_fnames:
        i_cat += 1
        aux_fname = '%s/%s' % (path_files, name)
        ccd = ccds[i_cat]
        mag_data = np.loadtxt(aux_fname, dtype=input_dt, usecols=input_icols, delimiter=catalog_delimiter)
        log.info('Reading file %s...' % aux_fname)
        if extinction:
            extinction_data = np.loadtxt(aux_fname, dtype=extinction_dt, usecols=extinction_icols,
                                         delimiter=catalog_delimiter)
            for key in extinction.keys():
                mag_data[key] = mag_data[key] - extinction_data[key]
        prop_data = np.loadtxt(aux_fname, dtype=properties_dt, usecols=properties_icols, delimiter=catalog_delimiter)

        try:
            Ngal = len(mag_data)
        except TypeError:  # It breaks if len == 1 :(
            Ngal = 1
            mag_data = [mag_data]

        ###
        ###  2 - Save redshift information (if exists)
        ###
        if z_col:
            redshift_data = np.loadtxt(aux_fname, dtype=np.float, usecols=(z_col,), delimiter=catalog_delimiter)
            db_redshift = db.get('/%s/%s/tables/z' % (filterset_name, ccd))
            if db_redshift is None:
                db_redshift = db.create_dataset(name='/%s/%s/tables/z' % (filterset_name, ccd),
                                                data=redshift_data,
                                                dtype=np.float,
                                                chunks=True,
                                                maxshape=(None,))
            else:
                db_redshift.resize((NgalCum[ccd] + Ngal,))
                db_redshift[NgalCum[ccd]:NgalCum[ccd] + Ngal] = redshift_data

        ###
        ###  3 - Save properties table data.
        ###
        db_properties = db.get('/%s/%s/tables/properties' % (filterset_name, ccd))
        if db_properties is None:
            db_properties = db.create_dataset(
                name='/%s/%s/tables/properties' % ((filterset_name, ccd)),
                data=prop_data,
                dtype=properties_dt,
                chunks=True,
                maxshape=(None,))
        else:
            db_properties.resize((NgalCum[ccd] + Ngal,))
            db_properties[NgalCum[ccd]:NgalCum[ccd] + Ngal] = prop_data

        ###
        ###  4 - Load magnitudes to the table.
        ###
        db_data = db.get('/%s/%s/data' % (filterset_name, ccd))
        if db_data is None:
            db_data = db.create_dataset(name='/%s/%s/data' % ((filterset_name, ccd)),
                                        shape=(Ngal, Nfilters),
                                        dtype=np.dtype([('m_ab', np.float), ('e_ab', np.float)]),
                                        chunks=(1, 5),
                                        maxshape=(None, Nfilters))
            NgalCum[ccd] = 0
        else:
            db_data.resize(((NgalCum[ccd] + Ngal, Nfilters)))

        for i_obj in range(Ngal):
            x = np.zeros(Nfilters, dtype=np.dtype([('m_ab', '<f8'), ('e_ab', '<f8')]))
            x['m_ab'] = [mag_data[i_obj][fil] for fil in filter_ids]
            x['e_ab'] = [mag_data[i_obj]['d%s' % fil] for fil in filter_ids]
            badpxl = np.bitwise_or(x['m_ab'] >= 97, x['m_ab'] <= -97)  #FIXME: Should be passed by the config file.
            x['m_ab'][badpxl] = np.inf
            x['e_ab'][badpxl] = np.inf

            db_data[i_obj + NgalCum[ccd]] = x

        NgalCum[ccd] += Ngal

    db.close()
