'''
Created on Sep 18, 2012

@author: william
'''

import numpy as np

from ..core.exceptions import MAGALException

#TODO: Remove this to use numpy.lib.recfunctions instead.
# http://projects.scipy.org/numpy/browser/trunk/numpy/lib/recfunctions.py?rev=8032 

def matchobjs(list1, list2):
    '''
        Create an id_list array with elements in list2 that matches list1.
        Useful to join tables by unique ids.
    '''
    id_list = []
    out_array = []
    list2 = np.sort(list2)
    aux_size2 = len(list2)
    for i_list in range(len(list1)):
        match = np.searchsorted(list2, list1[i_list])
        if (list1[i_list] == list2[match] and match < aux_size2):
            id_list.append(match)
            out_array.append(list2[match])
        else:
            raise MAGALException('Error. There are elements on list2 that does not exists on list1.')
            
    return id_list 

def get_zslice(l, z):
    i_z = np.argmin((l.z - z)**2) # Due to precision problems!
    return np.copy(l.library[i_z,:])