'''
Created on Feb 23, 2012

@author: william
'''

import os
import numpy as np
import h5py


class FilterSet(object):
    """
    This class reads a filterset file from file and returns a "filter" object.
    """

    def __init__(self, filterfile):
        '''
        Initialize a FilterSet object given a filterfile.

        Parameters
        ----------
        filterfile : string
            Path to the filter file.

        '''
        if not os.path.exists(filterfile):
            raise Exception('File not found: %s' % filterfile)

        if filterfile.endswith('.hdf5'):
            self._db_f = h5py.File(filterfile, 'r')
            self.filtersets = {}
            for f_ in self._db_f.keys():
                self.filtersets[f_] = [ccd for ccd in self._db_f[f_]]
        else:
            raise Exception('Unsupported filterfile.')


    def load(self, filterset_id, ccd):
        """
        Loads into ``self.filterset`` the ``filterset_id``
        
        Parameters
        ----------
        filterset_id : string
            Filterset identificator on filterfile. e.g.: ``sdss``

        ccd : string
            Identificator of which CCD should we load.
        
        Examples
        --------
                   
        See Also
        --------
        
        Notes
        -----
        """
        aux = self._db_f.get('/%s/%s' % (filterset_id, ccd))
        shape = sum([aux[m].len() for m in aux.keys()])
        self.filterset = np.empty(shape, dtype=np.dtype([('ID_filter', 'S32'), ('wl', np.float), ('transm', np.float)]))
        i = 0
        for id_filter in aux.keys():
            i_last = i + aux[id_filter].len()
            self.filterset['ID_filter'][i:i_last] = id_filter
            self.filterset['wl'][i:i_last] = aux[id_filter]['wl']
            self.filterset['transm'][i:i_last] = aux[id_filter]['transm']
            i = i_last

    @property
    def filterset_uniform(self, dl=1):
        """
        Modifies filter curves to match a specific uniform lambda coverage.

        Parameters
        ----------
        dl: float, optional
            Delta lambda spacing in Angstroms. (Default: 1 :math`\AA`)

        Returns
        -------
        filterset : array
            Interpolated filterset. Haves the same shape and dtype of ``FilterSet.filterset``
        """
        aux = []
        for fid in np.unique(self.filterset['ID_filter']):
            xx = self.filterset[self.filterset['ID_filter'] == fid]
            new_lambda = np.arange(xx['wl'].min(), xx['wl'].max(), 1.0)
            new_transm = np.interp(new_lambda, xx['wl'], xx['transm'])
            for i in range(len(new_lambda)):
                aux.append((fid, new_lambda[i], new_transm[i]))
        return np.array(aux, dtype=self.filterset.dtype)

    @property
    def filter_wls(self):
        '''
        Calculate the mean wave length of each filter. Useful for plotting.
        '''
        try:
            aux_names = np.unique(self.filterset['ID_filter'])
            dt = np.dtype([('ID_filter', 'S32'), ('wl_central', np.float)])
            aux_ret = np.array(
                [(fid, np.average(self.filterset[self.filterset['ID_filter'] == fid]['wl'])) for fid in aux_names],
                dtype=dt)
        except AttributeError:
            raise Exception('You have to load the filtersystem/CCD first!')

        return np.array(aux_ret)