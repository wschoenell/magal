import h5py

import numpy as np

from ..core.exceptions import MAGALException

from ..core.log import logger


class Library(object):
    '''
    Reads a .hdf5 template Library. 
    '''

    def _lib_open(self, library_file):
        # Open HDF5 file
        try:
            self._lib = h5py.File(library_file, 'r')
        except IOError:
            raise MAGALException('Could not open file %s.' % library_file)
        self._lib_filename = library_file

    def __init__(self, library_file):
        """
        Reads a .hdf5 template Library.
        
        Parameters
        ----------
        library_file: string
                      Template library file.
        """

        self.log = logger(__name__)

        self._lib_open(library_file)

        self.filtersystems = []
        self.ccds = {}
        for fid in self._lib.keys():
            if fid not in ('tables', 'ini_file'):
                self.filtersystems.append(fid)
                self.ccds[fid] = self._lib['/%s' % fid].keys()
        self.z = np.copy(self._lib['/tables/z'].value)
        self.properties = np.copy(self._lib['/tables/properties'].value)
        self._lib_isclosed = False

    def get_filtersys(self, fsys, ccd=None, close=True):
        """
        Select a filter system and a ccd on the library.
        Will raise an exception if ccd or filtersystem id does not exists.
        
        Parameters
        ----------
        fsys: string
              Filter system id.
              
        ccd: string
             CCD id.
        """
        if ccd is None:
            ccd = self.ccds[fsys][0]
        if fsys in self.filtersystems:
            if self._lib_isclosed:
                self._lib_open(self._lib_filename)
            self.fsys = fsys
            self.filterset = np.copy(self._lib['/%s/%s/filterset' % (fsys, ccd)].value)
            self.filtercurves = np.copy(self._lib['/%s/%s/filtercurves' % (fsys, ccd)].value)
            self.library = np.copy(self._lib['/%s/%s/library' % (fsys, ccd)].value)
            self.path = '/%s/%s/' % (fsys, ccd)
            self.log.debug('Read filtersystem %s/%s' % (fsys, ccd))
            self._lib.close()
            self._lib_isclosed = True
        else:
            raise MAGALException('Filtersystem/CCD %s/%s does not exists!' % (fsys, ccd))