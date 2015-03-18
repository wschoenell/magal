from magal.core.exceptions import MAGALException

__author__ = 'william'

from magal.library.ssp import SSP
from magal.library.csp import TwoExponential
from magal.library.starlight import StarlightSDSS

class LibraryModel(object):
    '''
    Abstract class for Template Library models manipulation.
    Used by magal_library script to create different types of libraries.
    '''

    def __init__(self, type, *kwargs):
        if type.startswith('two_exp'):
            self.__class__ = TwoExponential
            self.__init__(type, *kwargs)
        elif type == 'starlight_sdss':
            self.__class__ = StarlightSDSS
            self.__init__(type, *kwargs)
        elif type == 'ssp':
            self.__class__ = SSP
            self.__init__(type, *kwargs)
        else:
            raise MAGALException('Unknown library type %s' % type)