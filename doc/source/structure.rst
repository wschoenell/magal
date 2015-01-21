---------------------
MagAl Files Structure
---------------------

abc

*FilterSet* file
----------------

Contains transmission curves of the filters.

Attributes:
^^^^^^^^^^^

*None*

File structure:
^^^^^^^^^^^^^^^

``Filter_System`` is the filter system ID and ``CCD_1``, ``CCD_2``, etc are the CCD ids.

::

    .
    └── Filter_System
       ├── CCD_1
       |   ├── F_100
       |   ├── F_200
       |   └── ...
       └── CCD_2
           └── ...

*Input* file
------------

File containing the input data to ``magal_fit``. Can be either magnitudes, lick indices or both.

Attributes:
^^^^^^^^^^^

*None*

Template *Library* file
-----------------------

This file contains the models which will be used to fit the observables. Can be either magnitudes, lick indices or both.
This is used by ``magal_fit`` as template libraries and can also be used as input file when doing simulations.

Attributes:
^^^^^^^^^^^

*None*

File structure:
^^^^^^^^^^^^^^^

::

    .
    ├── Filter_System
    |   ├── CCD_1
    |   |   ├── filtercurves
    |   |   ├── filterset
    |   |   └── library
    |   └── CCD_2
    |       └── ...
    ├── tables
    |   ├── properties
    |   └── z
    └── ini_file


*Magal* Output file
-------------------

Attributes:
^^^^^^^^^^^

* ``input_file`` : Input filename
* ``input_md5`` : Input file MD5 sum
* ``library_file`` : Library filename
* ``library_md5`` : Library file MD5 sum

File structure:
^^^^^^^^^^^^^^^

``property_1``, ``property_1``, etc are the library properties from ``/tables/properties`` of ``Library`` file.
``Filter_System`` is the filter system ID and ``CCD_1``, ``CCD_2``, etc are the CCD ids from ``FilterSet`` file.

::

    .
    ├── Filter_System
    |   ├── CCD_1
    |   |   ├── chi2
    |   |   ├── likelihood
    |   |   ├── likelihood_template
    |   |   ├── n
    |   |   ├── s
    |   |   └── statistics
    |   |       └── statistics
    |   |           └── library
    |   |           |   ├── property_1
    |   |           |   |   ├── AVG
    |   |           |   |   └── percentiles
    |   |           |   ├── property_2
    |   |           |   |   ├── AVG
    |   |           |   |   └── percentiles
    |   |           |   └── ...
    |   |           └── model
    |   |               ├── i_BMX
    |   |               ├── n_percentiles
    |   |               └── s_Mass
    |   |                   ├── AVG
    |   |                   └── percentiles
    |   └── CCD_2
    |       └── ...
    └── ini_file
