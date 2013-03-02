'''
Created on Oct 22, 2012

@author: william
'''

import h5py

import logging
from magal.core.exceptions import MAGALException

class Library(object):
    '''
    Reads a .hdf5 template Library. 
    '''

    def __init__(self, library_file):
        '''
        Reads a .hdf5 template Library.
        
        Parameters
        ----------
        library_file: string
                      Template library file.
        '''
        
        self.log = logging.getLogger('magal.io.readinput')
        
        # Open HDF5 file
        try:
            self._lib =  h5py.File(library_file, 'r')
        except IOError:
            raise MAGALException('Could not open file %s.' % library_file)
        
        self.filtersystems = []
        self.ccds = {}
        for fid in self._lib.keys():
            if fid != 'tables':
                self.filtersystems.append(fid)
                self.ccds[fid] = self._lib['/%s' % fid].keys()
        self.z = self._lib['/tables/z'].value
        self.properties = self._lib['/tables/properties']
        
    def get_filtersys(self, fsys, ccd = None):
        '''
        Select a filter system and a ccd on the library.
        Will raise an exception if ccd or filtersystem id does not exists.
        
        Parameters
        ----------
        fsys: string
              Filter system id.
              
        ccd: string
             CCD id.
        '''
        if ccd == None: ccd = self.ccds[fsys][0]
        if fsys in self.filtersystems:
            self.filterset = self._lib['/%s/%s/filterset' % (fsys, ccd)].value
            self.filtercurves = self._lib['/%s/%s/filtercurves' % (fsys, ccd)].value
            self.library = self._lib['/%s/%s/library' % (fsys, ccd)]
            self.path = '/%s/%s/' % (fsys, ccd)
            self.log.debug('Read filtersystem %s/%s' % (fsys, ccd))
        else:
            raise MAGALException('Filtersystem/CCD %s/%s does not exists!' % (fsys, ccd))