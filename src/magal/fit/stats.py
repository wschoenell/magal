'''
Created on Oct 9, 2012

@author: william
'''

import numpy as np

def chi2(m_o, m_l, w):
    '''
    Returns n_good, s and chi2 for an observed-library AB magnitude pair.
    
    Parameters
    ----------        
              
    Returns
    -------
    n_good: array_like
            Number of good pixels
    
    s: array_like
       Scaling-factor. :math:`-2.5 \\log M_\\star = {\\sum w^2(l) * \\left(m_o(l) - m_l (l) \\right)}\\over{\\sum_l w^2(l)}`
       
    
    chi2: array_like
          Chi-square. :math:`\\chi^2 = \\sum_l \\left( m_o(l) - m_l(l) - s_{lo}(l) \\right)^2 * w^2(l)`
    
    Examples
    --------
               
    See Also
    --------
    
    Notes
    -----
    '''
    mask = (m_o == np.inf) + (m_l == np.inf) + (w == 0)
    mask = np.invert(mask) + np.isnan(m_o) + np.isnan(w)
    n_good = np.sum(mask)
    w2 = np.power(w[mask],2)
    s = np.sum( w2 * (m_o[mask] - m_l[mask]) ) / np.sum( w2 ) # Scaling-factor = - 2.5 log M_\star
    chi2 = np.sum( np.power(m_o[mask] - m_l[mask] - s,2) * w2 )
    
    return n_good, s, chi2

def percentiles(x,y,perc):
    '''
    Returns the qth percentile of the array elements to a given properties vector and values vetor.
    
    Parameters
    ----------
    x : array_like
        Properties values.

    y : array_like
        Array where percentiles are calculated on.
        
    
    perc : float in range of [0,100] (or sequence of floats)
           Percentile to compute which must be between 0 and 100 inclusive.
              
    Returns
    -------
    percentiles : float, ndarray #FIXME:
                  The interpolated values, same shape as `perc`.
                  
    See Also
    --------
    numpy.percentile
    
    '''
    perc = np.asanyarray(perc) / 100.
    aux_idx = np.argsort(x)
    y = y[aux_idx]
    x = x[aux_idx]
    y = np.cumsum(y)
    y = y/y[-1]
    out = np.interp(perc, y, x)
    return out