Fitting
=======

To run MagAl, you will need three files:

* Library File, generated by ``magal_library` script. See `libraries` for more info.

* Input File, which will be generated by the ``magal_input`` script. See `input` for more info.

* Configuration File, which will be described on the next section.



Ini file keywords
^^^^^^^^^^^^^^^^^

Section : ``[FitGeneral]``

* ``input_file`` : string. Input filename.

* ``library_file`` : string. Library filename.

* ``filter_sys`` : string. Filter system name.

* ``ccd`` : int. CCD number id.

* ``output_file`` : string. Output filename

* ``Nobj_max`` : int. Optional. Max number of objects to fit. Useful for tests.

* ``allow_overwrite`` : bool. Optional. If ``True``, allows ``magal_fit`` to overwrite the ``output_file``.

* ``zp_error`` : float. Optional. Zero point error. It will be quadratically summed to the galaxies error.
  E.g.: :math:`e_i = \sqrt{e_i^2 + e_{zp}^2}` .

* ``allow_multiprocessing`` : bool. Default: ``False``. Optional. If ``True``, ``magal_fit`` will use multiprocessing.

* ``Nz`` :

* ``filter_include`` :

* ``filters_exclude`` :

* ``fudge_factor`` :

Section : ``[FitGeneral]`` . Optional.

* ``is_simulation`` : bool

* ``obj_z`` :

* ``mass_field`` :


Example: NONONON
^^^^^^^^^^^^^^^^

aa