import os
import numpy as np
from magal.core.exceptions import MAGALException


__author__ = 'william'


class Spectral(object):
    """
    A library defined by a set of spectrum files.
    """

    def __init__(self, input_file, input_path=None, input_columns=None):
        """
        A library defined by a set of spectrum files.

        Parameters
        ----------
        input_file : string
            File with the model spectra filenames and the properties of each one.

        input_path : string. Optional.
            Path to the ``input_file`` files.

        input_columns : dict. Optional.
            Dictionary with the names and the columns of the model properties. All the properties are treated as
            ``numpy.float`` .

        Returns
        -------

        Examples
        --------
        >>> Spectral(filename, {'age': 2, 'Z': 4, 'Ha': 16})

        """
        self.input_file = input_file
        self.input_path = input_path
        self.input_columns = input_columns
        self._check_input()

        # Open the input_file and get the input_data properties
        dt = np.dtype([('t0_young', np.float), ('tau_young', np.float), ('t0_old', np.float), ('tau_old', np.float),
                       ('frac_young', np.float), ('tau_v', np.float), ('Z', np.float)])
        self.input_data = np.loadtxt(self.input_file, dtype=dt, usecols=(0,1,3))
        self.lib_size = len(self.input_data)

        pass

    def get_model_spectrum(self, i_model):
        """
        Return the spectrum of a given model on the models file.

        Parameters
        ----------
        i_model : int
            Model number which we want to take the spectrum.

        Returns
        -------
        spec : array
            An array with ``lambda`` and ``flux`` to the model ``i_model`` spectrum.
        """

        pass


def _check_input(self):
    # Sanity check.
    if self.type != 'spectral':
        raise MAGALException('Error creating Spectral LibraryModel object.')

    # Check if file exists
    if not os.path.isfile(self.input_file):
        raise IOError('File %s not found or is not a file.' % self.input_file)
    if not os.path.isdir(self.input_path):
        raise IOError('Directory %s not found or is not a directory.' % self.input_path)

    # Check if input_columns is a dict.
    if self.input_columns is not None and not isinstance(self.input_columns, dict):
        raise ValueError('input_columns is not a dictionary')
