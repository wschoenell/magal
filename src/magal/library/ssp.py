import os

import numpy as np
from pystarlight.util.base import StarlightBase

from magal.core.exceptions import MAGALException


class SSP(object):
    """
    Library from Starlight-like SSPs file
    """


    def __init__(self, library_type, base_file, base_path):
        """
        Library from Starlight-like SSPs file

        Parameters
        ----------
        library_type : string
            Will raise an error if it is different from ``two_exp``.

        base_file : string
            Starlight-like base filename.

        base_path : string
            Directory where the models files are. In case base_file is a hdf5 file, this is the path on hdf5 to the lib.

        """
        # 0 - Check if all the parameters are okay.
        self.type = library_type
        self.base_path = base_path
        self.base_file = base_file
        self._check_input()

        # 1 - Load data...
        # 1.1 - STARLIGHT base file + path.
        self.bt = StarlightBase(self.base_file, self.base_path)
        self.ages = np.unique(self.bt.ageBase)

        self.lib_size = self.bt.nAges * self.bt.nMet
        self._aux_metallicities, self._aux_ages = np.indices((self.bt.nMet, self.bt.nAges))
        self._aux_ages = self._aux_ages.ravel()
        self._aux_metallicities = self._aux_metallicities.ravel()

        # 1.2 - STARLIGHT base parameters
        dt = np.dtype([('age', np.float), ('Z', np.float), ('Mstars', np.float), ('aFe', np.float)])
        self.input_data = np.empty(self.lib_size, dtype=dt)
        for i_model in range(self.lib_size):
            self.input_data[i_model]['age'] = self.bt.ageBase[self._aux_ages[i_model]]
            self.input_data[i_model]['Z'] = self.bt.metBase[self._aux_metallicities[i_model]]
            self.input_data[i_model]['Mstars'] = self.bt.Mstars[self._aux_metallicities[i_model], self._aux_ages[i_model]]
            self.input_data[i_model]['aFe'] = self.bt.aFe[self._aux_metallicities[i_model], self._aux_ages[i_model]]

    def get_model_spectrum(self, i_model):
        """
        Return the spectrum of a given model on the models list.

        Parameters
        ----------
        i_model : int
            Model number which we want to take the spectrum.

        Returns
        -------
        spec : array
            An array with ``lambda`` and ``flux`` to the model ``i_model`` spectrum.
        """
        i_met, i_age = self._aux_metallicities[i_model], self._aux_ages[i_model]
        spec = np.zeros(shape=self.bt.f_ssp.shape[2], dtype=np.dtype([('wl', np.float), ('flux', np.float)]))
        spec['wl'] = self.bt.l_ssp
        spec['flux'] = self.bt.f_ssp[i_met, i_age]

        return spec[np.newaxis, :]  # For SSPs, like CSPs, we return only the model spectrum


    def _check_input(self):
        if self.type != 'ssp':
            raise MAGALException('Error creating SSP LibraryModel object.')

        aux_files = [self.base_file]

        for f_ in aux_files:
            if not os.path.isfile(f_):
                raise IOError('File %s not found or is not a file.' % f_)