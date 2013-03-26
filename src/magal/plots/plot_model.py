'''
Created on Feb 14, 2013

@author: william
'''

import numpy as np
from magal.fit.stats import percentiles as _percentiles

def plot_chi2(axis, prop, likelihood, x, x_correct=None, percentiles=None, nbins = 100, log=False):
    #plt.clf()
    
    if log:
        x = np.log10(x)
        x_correct = np.log10(x_correct)
        if percentiles:
            percentiles[i] = [np.log10(percentiles[i]) for i in percentiles.dtype.names]
    
    # Posterior #
    aux_hst1, aux_bins = np.histogram(x,weights=likelihood,bins=nbins, normed=True)
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts1 = left+(right-left)/2
    aux_hst1 = aux_hst1 / np.max(aux_hst1)
    
    axis.plot(pts1,aux_hst1, color='magenta')
    
    # Prior #
    aux_hst2, aux_bins = np.histogram(x, bins=nbins, normed=True)
    left = np.array(aux_bins[:-1])
    right = np.array(aux_bins[1:])
    pts2 = left+(right-left)/2
    aux_hst2 = aux_hst2 / np.max(aux_hst2)
    
    axis.plot(pts2,aux_hst2, ':', color='blue')
    
    ylim = (0., 1.01)
    if x_correct:
        axis.plot([x_correct,x_correct],ylim, color='green')
    axis.set_ylim(ylim)
    
    # Percentiles:
    if percentiles:
        if log:
            perc = [2.5, 16, 50, 84, 97.5]
            perc_dt = np.dtype([ ('%s' % p, np.float) for p in perc ])
            percentiles = np.array(tuple(_percentiles(x, likelihood, perc)), dtype=perc_dt)
        axis.plot([percentiles['16'], percentiles['16']],ylim, '-.', color='black')   #P16
        axis.plot([percentiles['84'], percentiles['84']],ylim, '-.', color='black')   #P84
        axis.plot([percentiles['50'], percentiles['50']],ylim, '--', color='red')   #P50
        
    
    xtxt = axis.get_xlim()[0] + (axis.get_xlim()[1] - axis.get_xlim()[0]) *.03
    ytxt = .72
    
    axis.text(xtxt, ytxt, prop)
    if x_correct:
        axis.plot([x_correct,x_correct],ylim, color='green')   #"Correct" value - To use with simulations
        axis.text(xtxt, ytxt-.1, '%5.3f' % (x_correct - percentiles['50']))


#def plot_model(magal_out):
