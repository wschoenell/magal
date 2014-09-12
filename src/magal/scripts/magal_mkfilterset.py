"""
Copyright (c) 2014, W. Schoenell
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the MagAl Team nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL W. SCHOENELL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
import sys
import time
import h5py
import atpy

import numpy as np

from magal.io.readfilterset import FilterSet
from magal.io.hdf5util import inithdf5

# if __name__ == '__main__' and len(sys.argv) > 2:

dbfile = sys.argv[1]
# Init file
db = inithdf5(dbfile)


for filter_file in sys.argv[2:]:
    aux_id = os.path.basename(filter_file).split('.')
    f = FilterSet()
    f.read(filter_file)
    for fid in np.unique(f.filterset['ID_filter']):
        dataset = '/%s/%s/%s' % (aux_id[0], aux_id[1] ,fid)
        print dataset
        aux = atpy.Table(name = fid)
        aux.add_column(name='wl', data = f.filterset['wl'][f.filterset['ID_filter'] == fid])
        aux.add_column(name='transm', data = f.filterset['transm'][f.filterset['ID_filter'] == fid])
        db.create_dataset(dataset, data = aux.data)


db.close()

# else:
#     print 'Usage: %s filterdbfile.hdf5 FilterSet.CCD1.filter FilterSet.CCD2.filter ... FilterSet.CCD#.filter' % sys.argv[0]