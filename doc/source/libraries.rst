----------------------
Making Magal Libraries
----------------------

bla bla bla

* Run ``magal_library library.ini``

Ini file keywords
-----------------

Section : ``[LibraryGeneral]``

* ``library_file`` : string. Output filename.
* ``filter_file`` : string. MagAl Filter transmission curves input file.
* ``filterset`` : string. Specifies which filter inside the filter_file will
  you use to create the library.
* ``library_type`` : string. Library type. Allowed types are: ``spectral``,
  ``two_exp`` and ``starlight_sdss`` .
* ``z_from`` : float. Minimum redshift on the library.
* ``z_to`` : float. Maximum redshift on the library.
* ``z_step`` : float. Delta-z of the library redshift binning.
* ``z_error`` : float. Optional. Random error :math:`\sigma` which will be summed to the spectroscopic redshift:
  :math:`z_{\rm photo} = z_{\rm spec} + \sigma (1 + z_{\rm spec})`
* ``allow_overwrite`` : bool. Optional. Default: False. If true, allow ovewriting output
  files.
* ``comsmology`` : dict. Optional. Default: WMAP9 (Hinshaw et al. 2013). Keys are ``H0``
  (:math:`H_0`) and ``Om0`` (:math:`\Omega_M`). Example:
  ``cosmology = {'H0': 70, 'Om0': 0.3}``


Simple library: Library from a filelist of spectra
--------------------------------------------------

On ``[LibraryGeneral]`` section, this library haves ``library_type = spectral``

Ini file keywords
^^^^^^^^^^^^^^^^^

Section : ``[LibraryParameters]``

* ``input_file`` : string. Input filename. Contains on the first column the ``filename`` of the
  spectrum and on the other columns the library properties.
* ``input_path`` : string. Optional. Input files path. Path to the be prepended to the ``filename`` .
* ``input_columns`` : dictionary. Optional. Contains the properties of each spectrum. Each property
  will be treated as ``numpy.float`` . Example: ``{'age': 2, 'Z': 6, 'log_Ha': 14}`` .

Example
^^^^^^^

Download: :download:`test_lib_spectral.ini <ini/test_lib_spectral.ini>`

.. literalinclude:: ini/test_lib_spectral.ini
   :language: ini


Parametric Library: Two-exponential parametric library
------------------------------------------------------

On ``[LibraryGeneral]`` section, this library haves ``library_type = two_exp``

Ini file keywords
^^^^^^^^^^^^^^^^^

Section : ``[LibraryParameters]``

* ``input_file`` : string. Input filename. Should contain 7 ``float``-type columns:
  ``t0_young, tau_young, t0_old, tau_old, frac_young, tau_v, Z``
* ``bases_dir`` : string. Directory prefix where the base model files are stored.
* ``base_file`` : string. STARLIGHT basefile location. For more details on base_file syntax
  consult the STARLIGHT Manual: http://starlight.ufsc.br/
* ``fraction_type`` : string. Type of weighting the model components: Allowed types ``light``
  and ``mass``.
* ``lambda_norm`` : float. Default: 4020 :math:`\unicode{x212B}`

Example
^^^^^^^

Download: :download:`test_lib_twoComponent.ini <ini/test_lib_twoComponent.ini>`

.. literalinclude:: ini/test_lib_twoComponent.ini
   :language: ini


STARLIGHT Library: STARLIGHT-SDSS project library
-------------------------------------------------

On ``[LibraryGeneral]`` section, this library haves ``library_type = starlight_sdss``

Ini file keywords
^^^^^^^^^^^^^^^^^

Section : ``[LibraryParameters]``

* ``input_file`` : string. Input filename. Should contain 2 ``string``-type columns:
  ``arq_obs, arq_syn``. For more details, consult STARLIGHT and its database documentation.
* ``tables_dir`` : string. Directory prefix where the STARLIGHT table files are stored.
* ``input_dir`` : string. Directory prefix where the STARLIGHT input files are stored.
* ``output_dir`` : string. Directory prefix where the STARLIGHT output files are stored.

Example
^^^^^^^

Download: :download:`test_lib_starlight.ini <ini/test_lib_starlight.ini>`

.. literalinclude:: ini/test_lib_starlight.ini
   :language: ini
