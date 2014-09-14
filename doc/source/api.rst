===
API
===


Contents:

.. contents::
   :local:
   :depth: 1


:mod:`magal.io`
===============

Classes used mainly to read MagAl output files.

.. automodule:: magal.io.readfilterset
   :members:

.. automodule:: magal.io.readinput
   :members:

.. automodule:: magal.io.readlibrary
   :members:

.. automodule:: magal.io.readmagal
   :members:



:mod:`magal.library`
====================

Classes used to make different kinds of libraries.

.. currentmodule:: magal.library.spectral
.. autoclass:: Spectral
   :members:

.. currentmodule:: magal.library.csp
.. autoclass:: TwoExponential
   :members:

.. currentmodule:: magal.library.starlight
.. autoclass:: StarlightSDSS
   :members:


:mod:`magal.fit`
================

.. currentmodule:: magal.fit.stats
.. autofunction:: chi2

.. currentmodule:: magal.fit.stats
.. autofunction:: percentiles


:mod:`magal.photometry.syntphot`
================================

.. currentmodule:: magal.photometry.syntphot
.. autofunction:: spec2filter

.. currentmodule:: magal.photometry.syntphot
.. autofunction:: spec2filterset

.. currentmodule:: magal.photometry.syntphot
.. autoclass:: photoconv
   :members:

:mod:`magal.plots`
==================

.. currentmodule:: magal.plots.plot_helpers
.. autoclass:: plotMagal
   :members:

.. currentmodule:: magal.plots.mosaic
.. autofunction:: get_mosaic

:mod:`magal.util`
=================

.. currentmodule:: magal.util.stellarpop
.. autoclass:: n_component
    :members:

.. currentmodule:: magal.util.cosmo
.. autofunction:: zcor

