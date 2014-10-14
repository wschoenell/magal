'''
Created on Jul 15, 2013

@author: william
'''

import numpy as np
import matplotlib.pyplot as plt
from magal.util.matchs import get_zslice
from ConfigParser import NoOptionError
import ast

def plot_magalFit(i, l, m, i_obj, ax=None, emission_lines=True, scale='linear'):
    
    
    # 3 - Plot fit model
    if ax == None:
        fig = plt.figure(1)
        fig.clf()
        ax = fig.add_axes((0.1,0.1,.8,.8))
            
    aux_l = get_zslice(l, i.z[i_obj])
    print i_obj

    # Check if only few filters were fitted.
    try:
        filters_include = ast.literal_eval(m.input_config.get('FitGeneral', 'filters_include'))
        filterset_mask = np.sum([k == l.filterset['ID_filter'] for k in filters_include], axis=0, dtype=np.bool)
    except NoOptionError:
        filters_include = None
    
    # Define few emission lines to plot.
    # TODO: This should be stored on other side of this shit.
    aux_labels = {'3727': '$[ O_{II} ]$', #  [\AA^{-1}]
             '4861': '$H\\beta$',
             '5007': '$O_{III}$',
             '6563': '$H_\\alpha$'} #,
#              '6584': '$\log W_{[ \\rm{N} II ]}$' }
    
    mask = (i.data.value[i_obj]['m_ab'] != np.inf)
    aux_model = np.tensordot(m.likelihood[i_obj], np.array([aux_l['m_ab'][i_prop] + m.s[i_obj][0][i_prop] for i_prop in range(len(aux_l))]), [1,0])[0]
    if scale != 'linear':
        ax.errorbar(l.filterset['wl_central'][mask], i.data.value[i_obj]['m_ab'][mask], yerr=i.data.value[i_obj]['e_ab'][mask], fmt=None, marker=None)
        ax.plot(l.filterset['wl_central'], aux_model, '.', color='red')
#         print i.data.value[i_obj]['m_ab'][mask]
        ax.plot(l.filterset['wl_central'][mask], i.data.value[i_obj]['m_ab'][mask], '.', color='blue')
        if filters_include:
            ax.plot(l.filterset['wl_central'][filterset_mask], aux_model[filterset_mask], '.', color='green')
    else:
        ax.errorbar(l.filterset['wl_central'][mask], 10**(-.4 * i.data.value[i_obj]['m_ab'][mask]), yerr=10**(-.4 * i.data.value[i_obj]['e_ab'][mask]) - 1, fmt=None, marker=None)
        ax.plot(l.filterset['wl_central'], 10**(-.4 * aux_model), '.', color='red')
#         print i.data.value[i_obj]['m_ab'][mask]
        ax.plot(l.filterset['wl_central'][mask], 10**(-.4 * i.data.value[i_obj]['m_ab'][mask]), '.', color='blue')
        if filters_include:
            ax.plot(l.filterset['wl_central'][filterset_mask], 10**(-.4 * aux_model[filterset_mask]), '.', color='green')

            

#     ax.plot(l.filterset['wl_central'], aux_model, '.', color='red')
#     aux_model = aux_l[m.i_BMX[i_obj]][0]
#     aux_model['m_ab'] += m.s[i_obj, m.i_BMX[i_obj][0], m.i_BMX[i_obj][1]] 
#     np.array([aux_l['m_ab'][i] + m.s[i_obj][0][i] for i in range(len(aux_l))]) 


    z_corr = (1.+0)/(1.+i.z[i_obj]) # (1.+toz)/(1.+fromz)
    aux_rest_lambda = ax.get_xticks() * z_corr
    aux_label = []
    for l_ in aux_rest_lambda:
        aux_label.append('%3.2f' % l_)
    ax_aux = ax.twiny()
    ax_aux.set_xlim(ax.get_xlim())
    ax_aux.set_xticklabels(aux_label)
    ylim = np.array(ax.get_ylim())
    ylim[0] = ylim[0]- (ylim[1] * .1) 
    ylim[1] = ylim[1]*1.1
    
    if emission_lines:
        for l_ in aux_labels.keys():
            ax.plot(np.divide([int(l_),int(l_)], z_corr), ylim, '--', color='black')
            ax.text(np.divide(int(l_), z_corr), ylim[1]*.9, aux_labels[l_], rotation='vertical', color='red')
        ax.set_ylim(ylim[0], ylim[1])