'''
Created on Oct 22, 2012

@author: william
'''

import h5py

import logging
from magal.core.exceptions import MAGALException

class Input(h5py.File):
    '''
    Reads a .hdf5 Input photometry file. 
    '''

    def __init__(self, photo_file):
        '''
        Reads a.hdf5 input photometry file.
        
        Parameters
        ----------
        photo_file: string
                    Photometry input file.
        '''
       

        log = logging.getLogger('magal.io.readinput')
        
        # 0 - Open HDF5 file
        h5py.File.__init__(self, photo_file, 'r')
        
        # 1.1 - Get Filtersets
        self.filtersystems = self.keys()
        self.filtersystems.remove('tables')
        # 1.2 - and CCDs
        self.ccds = [key for key in self[self.filtersystems[0]].keys()]
        
        try:
            self.z = self['/tables/properties'].value['z']
        except:
            log.warning('Warning: The input photometry file does not have redshift table!')
        self.properties = self['/tables/properties']
        
    def get_filtersys(self, fsys, ccd):
        '''
        Select a filter system and a ccd on the inputfile.
        Will raise an exception if ccd or filtersystem id does not exists.
        
        Parameters
        ----------
        fsys: string
              Filter system id.
              
        ccd: string
             CCD id.
        '''
        
        if fsys in self.filtersystems:
            self.data = self['/%s/%s/data' % (fsys, ccd)]
            self.path = '/%s/%s/' % (fsys, ccd)
        else:
            raise MAGALException('Filtersystem/CCD %s/%s does not exists!' % (fsys, ccd))