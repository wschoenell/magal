Making filtersets file
----------------------

To start playing around with **MagAl** you should have the filter curves stored on its format. To create it from ASCII
files, you should run the ``magal_mkfilterset`` script with the ``.ini`` configuration as described on this page.

Ini file keywords
^^^^^^^^^^^^^^^^^

Section : ``[FilterLibrary]``

* ``filterset_file`` : string. Output filename.
* ``filterset_id`` : string. Filterset identification.
* ``filter_names`` : dict. Dictionary where the key are the ccd names and
  the values are lists of filter names.
* ``filter_filenames`` : dict. Dictionary where the key are the ccd names
  and the values are lists of filter filenames. From the ascii filter files it will be 
  used the 2 first columns which must be wavelength (in :math:`\unicode{x212B}`)
  and transmission (in arbitrary units).
  
  
Example: Making filterset file to SDSS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


* Download the SDSS filter transmission curves and save, for example, on
  $HOME/Downloads/sdss/ folder. SDSS filters can be downloaded from this page:
  http://www.sdss2.org/dr7/instruments/imager/index.html#tables

* Make an .ini file for it:

Download: :download:`mkfilter.ini <ini/mkfilter.ini>`

.. literalinclude:: ini/mkfilter.ini
   :language: ini
   
* Run ``magal_mkfilterset mkfilter.ini``