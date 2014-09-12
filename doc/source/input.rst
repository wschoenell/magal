Making input file
-----------------

bla bla bla

Ini file keywords:

Section : ``[InputGeneral]``

* ``output_file`` : string. Output filename.
* ``filterset_id`` : string. Filterset name. Must be equal to the ``magal_mkfilter`` id.
* ``path_files`` : string. Path to the magnitude catalog files.
* ``cat_files`` : list. List of strings with the name of the catalog files.
* ``ccds`` : list.
* ``delimiter`` : string. Optional. 
* ``columns`` : dict. 
* ``extinction`` : dict or bool. Optional. If True, will calculate extinction with CCM89 reddening law and Schlegel dust maps.

Section : ``[InputProperties]``

* ``unique_id_column`` : int.
* ``general`` : dict.
* ``flags`` : dict. Optional.
* ``redshift`` : int. Optional.
  
  
Example: Making filterset file to SDSS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


* Download the SDSS filter transmission curves and save, for example, on
  $HOME/Downloads/sdss. SDSS filters can be downloaded from this page:
  http://www.sdss2.org/dr7/instruments/imager/index.html#tables

* Make an .ini file for it:

Download: :download:`mkfilter.ini <ini/mkfilter.ini>`

.. literalinclude:: ini/mkfilter.ini
   :language: ini
   
* Run ``magal_mkfilterset mkfilter.ini`` 