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
import sys

import numpy as np
from magal.io.hdf5util import inithdf5

from magal.core.config import MagalConfig


def main():
    if len(sys.argv) > 1:

        config = MagalConfig(sys.argv[1], 'mkfilterset')
        if not config:
            sys.exit(2)

        # Init file
        db = inithdf5(config.filterset_file)

        dt = np.dtype([('wl', np.float), ('transm', np.float)])

        filter_dset = None

        for ccd, filter_files in config.filter_filenames.iteritems():
            i_name = 0
            for filter_file in filter_files:
                filter_id = config.filter_names[ccd][i_name]
                dataset = '/%s/%s/%s' % (config.filterset_name, ccd, filter_id)
                filter_data = np.loadtxt(filter_file, dtype=dt)
                if filter_dset is None:
                    db.create_dataset(dataset, data=filter_data, dtype=dt)
                else:
                    print 'Error. Found duplicate filter_names.'
                    return 2
                i_name += 1

        db.close()
        return 0
    else:
        print 'Usage: %s configuration.ini' % sys.argv[0]
        return 2