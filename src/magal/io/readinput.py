import h5py

from ..core.exceptions import MAGALException
from ..core.log import logger


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
        
        # 0 - Open HDF5 file
        h5py.File.__init__(self, photo_file, 'r')
        
        # 1.1 - Get Filtersets
        self.filtersystems = self.keys()
        self.filtersystems.remove('ini_file')
        # 1.2 - and CCDs
        self.ccds = [key for key in self[self.filtersystems[0]].keys()]
        
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
        
        log = logger(__name__)
        
        if fsys in self.filtersystems:
            self.ccd = ccd
            self.fsys = fsys
            self.data = self['/%s/%s/data' % (fsys, ccd)]
            self.path = '/%s/%s/' % (fsys, ccd)
            try:
                self.z = self['/%s/%s/tables/z' % (fsys, ccd)].value
            except:
                log.warning('The input photometry file does not have redshift table!')
            try:
                self.properties = self['/%s/%s/tables/properties' % (fsys, ccd)].value
            except:
                log.warning('The input photometry file does not have properties table!')
        else:
            raise MAGALException('Filtersystem/CCD %s/%s does not exists!' % (fsys, ccd))
