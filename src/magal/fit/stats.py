'''
Created on Oct 9, 2012

@author: william
'''

import numpy as np


def chi2(m_o, m_l, w=None):
    """
    Returns n_good, s and chi2 for an observed-library AB magnitude pair.
    
    Parameters
    ----------
    m_o : array
        Observed magnitudes.

    m_l : array
        Library model magnitudes.

    w : array
        Chi-square weight. See below.
              
    Returns
    -------
    n_good : array_like
            Number of good pixels
    
    s : array_like
       Scaling-factor. :math:`-2.5 \\log M_\\star = {\\sum w^2(l) * \\left(m_o(l) - m_l (l) \\right)}\\over{\\sum_l w^2(l)}`
       
    
    chi2: array_like
          Chi-square. :math:`\\chi^2 = \\sum_l \\left( m_o(l) - m_l(l) - s_{lo}(l) \\right)^2 * w^2(l)`
    """
    if w is None:
        w = np.ones_like(m_o)
    mask = (m_o == np.inf) + (m_l == np.inf) + (w == 0) + np.isnan(m_o) + np.isnan(w)
    mask = np.invert(mask)
    w2 = w[mask]
    w2 **= 2
    dm = m_o[mask] - m_l[mask]
    s = np.sum(w2 * dm) / np.sum(w2)  # Scaling-factor = - 2.5 log M_\star

    n_total = len(mask)
    n_good = len(w2)
    if n_good == 0:
        print 'n_good is ZERO!'
    # return n_good, s, np.sum(np.power(dm - s, 2) * w2) # n_good, s, chi2
    return n_good, s, n_good * np.sum(np.power(dm - s, 2) * w2) / n_total  # n_good, s, chi2


def percentiles(x, y, perc):
    """
    Returns the qth weighted percentile of the array elements to a given properties vector and values vector.

    Parameters
    ----------
    x : array_like
        Properties values.

    y : array_like
        Array where percentiles are calculated on (weight).

    perc : float in range of [0,100] (or sequence of floats)
        Percentile to compute which must be between 0 and 100 inclusive.

    Returns
    -------
    percentiles : ndarray
        The interpolated values, same shape as ``perc``.

    See Also
    --------
    numpy.percentile
    """
    perc = np.asanyarray(perc) / 100.
    aux_idx = np.argsort(x)
    y = y[aux_idx]
    x = x[aux_idx]
    y = np.cumsum(y)
    y /= y[-1]
    return np.interp(perc, y, x)