.. _LES:

*************************************
Optimization with LES Data
*************************************

The main routine contained in ``runLESanalysis.py`` performs model optimization and bias correction using Large Eddy Simulation (LES) data from WRF-SFIRE. This routine is not required for running the plume rise model in predictive mode. Default bias parameters stored in ``config.py`` will be sourced for all general cases.

Key configuration setting for this analysis are contained in ``config.py``. Their descriptions can be found in :ref:`config-table`.


.. note::

   Running this part of the code requires access to interpolated cross-sectional data from WRF-SFIRE generated synthetic plumes.


Table of Configuration Settings
-------------------------------
.. _config-table:


.. list-table:: Configuration Settings
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
