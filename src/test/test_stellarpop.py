'''
Created on Dec 10, 2012

@author: william
'''

import atpy
import pystarlight.io
import matplotlib.pyplot as plt

from read_bc03bin import Read_BC03_Binary

from magal.util.stellarpop import n_component

from pystarlight.util.constants import L_sun

import numpy as np

#if __name__ == '__main__':
    
basedir = '/Users/william/mestrado/BasesDir/'
basefile = 'Base.bc03.Padova1994.chab.All'

bt = atpy.Table(basedir+basefile, basedir, read_basedir=True, type='starlightv4_base')

Z = np.unique(bt['Z_base']) 
ages = np.unique(bt['age_base'])

model = n_component(ages)
t_eb = 5.51000000e+07 #10**8.25
tau = 1e5 #10**8.0
frac = 1.0 # Fraction on exponential
model.add_exp(t_eb, tau, frac)
#    model.add_square(t0, l, 1-frac )

model.plot()

plt.figure(2)
plt.clf()

sfh = model.get_sfh()

norm = np.trapz(sfh, ages) #L[wl == 4020.0]

sfh = sfh.reshape(sfh.shape[0], 1)
ssps = bt.f_ssp[bt.Z_base == 0.02]
wl = bt.l_ssp[0]

L = np.trapz(sfh * ssps, ages, axis=0)

#Plot our model

aux_lim = np.bitwise_and(wl > 3000, wl < 9000)
plt.clf()
plt.plot(wl[aux_lim], L[aux_lim]/norm, color='red')

#Compare to an SSP with age = t_eb
norm = 1 #ssps[np.argwhere(ages < t_eb)[-1]][0][wl == 4020.0]
plt.plot(wl[aux_lim], ssps[np.argwhere(ages < t_eb)[-1]][0][aux_lim]/norm, color='cyan')


##Compare to BC03 code:
## FIXME: This does not work... I did not understand how are structured bc03 files. Must implement this someday.
#f_ised = 'data_example/bc03_stuff/csp_test.ised'
#bc03_data = Read_BC03_Binary(f_ised, big_endian=False, Full=True)
#bc03_wl = bc03_data[1]
#i_csp = int(np.argwhere(bc03_data[0] < 10**7.25)[-1])
#bc03_flux = bc03_data[2].T[i_csp]
#
#norm = L_sun #bc03_flux[bc03_wl == 4020.0]
#aux_lim = np.bitwise_and(bc03_wl > 3000, bc03_wl < 9000)
#plt.plot(bc03_wl[aux_lim], bc03_flux[aux_lim]/norm, color='green')

