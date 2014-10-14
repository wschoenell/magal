'''
Created on Jul 23, 2013

@author: william
'''

import numpy as np

from ..core.log import logger
from ..util.matchs import get_zslice2
from ..fit.stats import chi2


def chi2_wrapper(p):
    log = logger(__name__)
    i_obj, obj, filterset_mask, N_tpl, N_z, inp_z, is_simulation, Nz, lib, lib_z, mass_field, zp_error = p

    if i_obj % 10 == 0:
        log.debug('i_obj: %i' % i_obj)

    aux_n = np.empty(shape=(N_z, N_tpl))
    aux_s = np.empty(shape=(N_z, N_tpl))
    aux_chi2 = np.empty(shape=(N_z, N_tpl))
    for i_z in range(N_z):
        if Nz:  # If redshift comes from inputfile.
            a = get_zslice2(lib, lib_z, inp_z[i_obj])
        else:
            a = get_zslice2(lib, lib_z, inp_z[i_z])
        for i_tpl in range(N_tpl):
            # 3.2 - If this is a simulation, use 1% of magnitude as error.

            e_ab = obj['e_ab'][filterset_mask].copy()  #FIXME: ADD base error
            if zp_error:
                e_ab = np.sqrt(e_ab ** 2 + zp_error ** 2)
            w = 1 / e_ab
            del e_ab

            if is_simulation and mass_field is not None:
                aux_n[i_z, i_tpl], aux_s[i_z, i_tpl], aux_chi2[i_z, i_tpl] = chi2(
                    obj['m_ab'][filterset_mask] - mass_field[i_obj] * 2.5, a[i_tpl]['m_ab'][filterset_mask], w)
            else:
                aux_n[i_z, i_tpl], aux_s[i_z, i_tpl], aux_chi2[i_z, i_tpl] = chi2(obj['m_ab'][filterset_mask],
                                                                                  a[i_tpl]['m_ab'][filterset_mask], w)

    return i_obj, aux_n, aux_s, aux_chi2


def chi2_parameters(Nz, is_simulation, o_list, N_obj, N_z, N_tpl, inp_z, filterset_mask, lib, mass_field, zp_error):
    for i_obj in range(N_obj):
        obj = o_list[i_obj]
        yield i_obj, obj, filterset_mask, N_tpl, N_z, inp_z, is_simulation, Nz, lib.library, lib.z, mass_field, zp_error