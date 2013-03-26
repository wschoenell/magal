#!/bin/bash
# This script should be used to test magal main script functionalities...
db_dir="/Users/william/downloads/databases/"
out_dir="./"

rm test_real.hdf5 test_simu.hdf5 magal_real.profile magal_simul.profile

# Test 1: Real SDSS data.
python -m cProfile -o magal_real.profile ../scripts/magal -i $db_dir/sdss_modelMag.hdf5 -l $db_dir/lib_csp_A_sdss_z0p7.hdf5 -f sdss -c 1 -o test_real.hdf5 -N 3 -nz -v

# Test 2: Simulate with a Library as input.
python -m cProfile -o magal_simul.profile ../scripts/magal -i $db_dir/lib_csp_A_sdss_z0p7.hdf5 -l $db_dir/lib_csp_A_sdss_z0p7.hdf5 -f sdss -c 1 -o test_simu.hdf5 -N 1 -s -z 0.1 -v
