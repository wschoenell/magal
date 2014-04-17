'''
Created on Aug 13, 2012

@author: william
'''

import numpy as np

def zcor(spec, toz, fromz=0.0):
    '''
    Shift a spectrum from one to another redshift.
        
    Parameters
    ----------
    spec: dict
          Spectrum
          { wl: Wavelength (in Angstroms!)
            flux: Flux on a given wl
            error: Flux error. (Optional) }
    toz: float
         Redshift in which we want the output spectrum
    fromz: float
           CURRENT redshift of the spectra (default: 0.0, rest-frame)
    
    Returns
    -------
    s: dict
       Shifted spectrum
    '''
    s = np.copy(spec)
    k = (1.+toz)/(1.+fromz)
    s['wl'] = s['wl'] * k
    k = 1./k
    s['flux'] = s['flux'] * k
    if('error' in s.dtype.names):
        s['error'] = spec['error'] * k
    return s