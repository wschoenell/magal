Class Example:
--------------

Abstract class for Template Library models manipulation.
Used by magal_library script to create different types of libraries.

See also
--------
fucntion1, function2


Function Example:
-----------------


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