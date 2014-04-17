'''
Created on Jul 23, 2013

@author: william
'''

import numpy as np
from magal.util.matchs import get_zslice
from magal.fit.stats import chi2
from magal.io.readlibrary import Library

def chi2_wrapper(p):
#     i_obj,i_z,i_tpl, o, m, w = p
    i_obj, obj, filterset_mask, N_tpl, N_z, inp_z, is_simulation, Nz, lib_file, filter_sys, ccd, mass_field = p
    
    if i_obj % 10 == 0:
        print 'i_obj: %i' % i_obj
    
    # 1.2 - Library
    lib = Library(lib_file) #config.get('FitGeneral', 'lib_file'))
    lib.get_filtersys(filter_sys, ccd) #config.get('FitGeneral', 'lib_file')), config.get('FitGeneral', 'filter_sys'), config.get('FitGeneral', 'ccd')) # libra

    aux_n = np.empty(shape=(N_z,N_tpl))
    aux_s = np.empty(shape=(N_z,N_tpl))        
    aux_chi2 = np.empty(shape=(N_z,N_tpl))
    for i_z in range(N_z):
        if Nz: # If redshift comes from inputfile.
            a = get_zslice(lib, inp_z[i_obj])
#             log.debug('Redshift slicing inp_z, lib_z: %3.4f, %3.4f' % (inp_z[i_obj], lib.z[np.argmin((lib.z - inp_z[i_obj]) ** 2)]))
        else:
            a = get_zslice(lib, inp_z[i_z])
        for i_tpl in range(N_tpl):
        # 3.2 - If this is a simulation, use 1% of magnitude as error.
            if (is_simulation):
                w = 1 / (obj['e_ab'][filterset_mask] + 0.05) #TODO: Put a more realistic error component.
            else:
                w = 1 / obj['e_ab'][filterset_mask]

            if is_simulation and mass_field is not None:
                aux_n[i_z, i_tpl], aux_s[i_z, i_tpl], aux_chi2[i_z, i_tpl] = chi2(obj['m_ab'][filterset_mask] - 2.5 * mass_field[i_obj], a[i_tpl]['m_ab'][filterset_mask], w)
            else:
                aux_n[i_z, i_tpl], aux_s[i_z, i_tpl], aux_chi2[i_z, i_tpl] = chi2(obj['m_ab'][filterset_mask], a[i_tpl]['m_ab'][filterset_mask], w)

    return i_obj, aux_n, aux_s, aux_chi2

def chi2_parameters(Nz, is_simulation, o_list, N_obj, N_z, N_tpl, inp_z, filterset_mask, lib_file, filter_sys, ccd, mass_field, e_a_veia_a_fiar=None):
    for i_obj in range(N_obj):
        obj = o_list[i_obj]
        yield i_obj, obj, filterset_mask, N_tpl, N_z, inp_z, is_simulation, Nz, lib_file, filter_sys, ccd, mass_field

    
#                 n_ds[i_obj,i_z,i_tpl], s_ds[i_obj,i_z,i_tpl], chi2_ds[i_obj,i_z,i_tpl] = chi2(obj['m_ab'][filterset_mask], a[i_tpl]['m_ab'][filterset_mask], w)
#                 if (i_z % 5 == 0 or config.getboolean('FitGeneral', 'Nz')) and i_tpl == 0 and args.verbose > 0:
#                     log.setLevel('DEBUG')
#                     log.debug('I\'m at i_obj, i_z --> %s, %s' % (i_obj, i_z))
#                 if args.verbose > 0 and i_tpl == 0:
#                     print 'I\'m at i_obj, i_z --> %s, %s' % (i_obj, i_z)