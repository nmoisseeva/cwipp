How To Use
============

Introduction
---------------

The contents of this package include two process-based routines, several methods implementing CWIPP algorithms and a collection of plotting functions.

The key routines are summarized below:

- ``runCWIPP.py`` - Runs the plume rise model for a general case (in predictive mode). Steps are  summarized in :ref:`basic-usage` below.

- ``runLESanalysis.py`` - Performs optimization of CWIPP model using synthetic WRF-SFIRE LES plume data. Details can be found :ref:`here <LES>` (this routine is optional).



.. _basic-usage:

Basic Usage
--------------


``runCWIPP.py`` contains sample code for applying the CWIPP model in a forecast setting to predict the vertical smoke profile of a real-world wildfire plume

Key configuration setting for this analysis are contained in ``config.py``. Their descriptions can be found in :ref:`config-table`.

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

Configuration Settings
-------------------------------
.. _config-table:


.. list-table::
   :widths: 30 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - **zstep**
     - int
     - height interpolation step for analysis
   * - **interpZ**
     - 1D array
     - vector corresponding to vertical levels for analysis
   * - **BLfrac**
     - float
     - fraction of BL height to use as reference height z_s (default = 0.75)
   * - **g**
     - float
     - gravity constant (9.81 m/s)
   * - **PMcutoff**
     - float
     - minimum PM value to define plume edge
   * - **biasFit**
     - list
     - default model bias fit parameters ([0.9195, 137.9193])
   * - **figdir**
     - str
     - Set path for storing figures
