How To Use
============

Introduction
---------------

The contents of this package include two variable computation/extraction
routines, several methods implementing CWIPP algorithms and a collection of plotting functions.

The key routines are summarized below:

- ``runCWIPP.py`` - Runs the plume rise model for a general case (in predictive mode). General steps are  summarized in :ref:`basic-usage` below.

- ``runLESanalysis.py`` - Performs optimization of CWIPP model using synthetic WRF-SFIRE LES plume data. Details can be found :ref:`here <LES>`.



.. _basic-usage:

Basic Usage
--------------


``runCWIPP.py`` contains sample code for applying the CWIPP model in a forecast setting to predict the vertical smoke profile of a real-world wildfire plume

.. .. code-block:: python
..
..     from __future__ import print_function
..
..     from netCDF4 import Dataset
..     from wrf import getvar
..
..     ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
..
..     # Get the Sea Level Pressure
..     slp = getvar(ncfile, "slp")
..
..     print(slp)

Result:
