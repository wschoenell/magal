rm -i test_sdss.hdf5
python ../scripts/sdss2photdb.py -i ./data_example/test_sdss.txt -o test_sdss.hdf5 -f ./data_example/JPAS_51.hdf5 -v
