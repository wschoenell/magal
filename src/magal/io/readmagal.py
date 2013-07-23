'''
Created on Feb 14, 2013

@author: william
'''

import h5py

from magal.core.exceptions import MAGALException
import os
import ConfigParser

class MagalFit(h5py.File):
    '''
    Reads a .hdf5 MagalFit file. 
    '''

    def __init__(self, magal_file):
        '''
        Reads a MagalFit output file.
        
        Parameters
        ----------
        photo_file: string
                    Photometry input file.
        '''
        
        # 0 - Open HDF5 file
        try:
            h5py.File.__init__(self, magal_file, 'r')
        except IOError:
            raise MAGALException('Could not open file %s.' % magal_file)
        
        # 1.1 - Get Filtersets
        self.filtersystems = self.keys()
        
        # 1.2 - Get CCDs
        self.ccds = [key for key in self[self.filtersystems[0]].keys()]
        
        # 1.3 - Get stat properties
        self.props = self['%s/%s/statistics/library/' % (self.filtersystems[0], self.ccds[0])].keys() 
        self.stats = self['%s/%s/statistics/library/%s' % (self.filtersystems[0], self.ccds[0], self.props[0])].keys() 
        
        # 1.4 - Def some useful data
        self.N_prop = len(self.props)
        
        # 1.5 - Read configuration if exists
        try:
            self.ini_file = self['/ini_file'][0]
            a = os.tmpfile()
            a.write(self.ini_file)
            a.seek(0)
            self.input_config = ConfigParser.ConfigParser()
            self.input_config.readfp(a)
            a.close()
        except KeyError: # If configuration does not exists, move on...
            self.ini_file = ''
            self.input_config = ConfigParser.ConfigParser()
        
        
    
    def get_filtersys(self, fsys, ccd):
        '''
        Select a filter system and a ccd on the magal fitfile.
        Will raise an exception if ccd or filtersystem id does not exists.
        
        Parameters
        ----------
        fsys: string
              Filter system id.
              
        ccd: string
             CCD id.
        '''
        
        if fsys in self.filtersystems:
            
            # Fit data:
            self.chi2 = self['/%s/%s/chi2' % (fsys, ccd)]
            self.n = self['/%s/%s/n' % (fsys, ccd)]
            self.s = self['/%s/%s/s' % (fsys, ccd)]    
            
            # Useful data: TODO: READ THIS FROM FILE attrs!
            self.N_obj = self.chi2.shape[0]
            self.N_z = self.chi2.shape[1]  
            
            self.likelihood = self['/%s/%s/likelihood' % (fsys, ccd)]
            self.likelihood_T = self['/%s/%s/likelihood_template' % (fsys, ccd)]
            if self.N_z > 1:
                self.likelihood_z = self['/%s/%s/likelihood_redshift' % (fsys, ccd)]
            
            # Library statistics:
            for prop in self.props:
                self.__setattr__(prop, self['/%s/%s/statistics/library/%s' % (fsys, ccd, prop)]) # Properties
                for stat in self.stats:
                    self.__getattribute__(prop).__setattr__(stat, self['/%s/%s/statistics/library/%s/%s' % (fsys, ccd, prop, stat)]) # Their statistics
            
            # Model statistics
            self.i_BMX = self['/%s/%s/statistics/model/i_BMX' % (fsys, ccd)]
            self.s_Mass = self['/%s/%s/statistics/model/s_Mass' % (fsys, ccd)]
            self.s_Mass.AVG = self['/%s/%s/statistics/model/s_Mass/AVG' % (fsys, ccd)]
            self.s_Mass.percentiles = self['/%s/%s/statistics/model/s_Mass/percentiles' % (fsys, ccd)]
            
            self.path = '/%s/%s/' % (fsys, ccd)
        else:
            raise MAGALException('Filtersystem/CCD %s/%s does not exists!' % (fsys, ccd))