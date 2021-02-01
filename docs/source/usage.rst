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

``runCWIPP.py`` contains sample code for applying the CWIPP model in a forecast setting to predict the vertical smoke profile of a real-world wildfire plume.

Key configuration setting for this analysis are contained in ``config.py``. Their descriptions are summarized in Configuration Settings below.


Input Data
--------------
Sample data for testing can be found in ``sample_fires.npy``.

The model expects dictionary input, where each case is identified by plume name/ID (string). Each plume in the parent input dictionary is also a dictionary, containing information fireline intensity ``I`` and vertical potential temperature profile ``sounding``.

.. code-block:: python

  >>> inputs[case.name]
  {'sounding': array([...]),
  'I': float,
  'zCL': float}

After successful completion of ``runCWIPP.py`` the same dictionary is appended with two additional keys (modelled mean injection height ``zCL`` and classification ``penetrative``) and saved as ``output.npy``

.. code-block:: python

  >>> inputs[case.name]
  {'sounding': array([...]),
  'I': float,
  'zCL': float,
  'zCL': float,
  'penetrative': boolean}


Configuration Settings
-------------------------------
.. _config-table:


.. list-table::
   :widths: 30 10 60
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - **input_data**
     - dict
     - input information about fires (id, soundings, intensity)
   * - **dz**
     - int
     - vertical step (in meters) between sounding levels
   * - **zstep**
     - int
     - height interpolation step for analysis (in meters)
   * - **interpZ**
     - 1D array
     - vector corresponding to vertical levels for analysis
   * - **BLfrac**
     - float
     - fraction of BL height to use as reference height z_s (default = 0.75)
   * - **biasFit**
     - list
     - default model bias fit parameters ([0.9195, 137.9193])
   * - **figdir**
     - str
     - Set path for storing figures
