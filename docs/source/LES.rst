.. _LES:

*************************************
CWIPP Optimization with LES Data
*************************************

The main routine contained in ``runLESanalysis.py`` performs model optimization and bias correction using Large Eddy Simulation (LES) data from WRF-SFIRE. This routine is not required for running the plume rise model in a predictive mode. Default bias parameters stored in ``config.py`` will be sourced for all general cases. 

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
   * - **wrfdir**
     - str
     - Path to interpolated cross-sectional plume data
   * - **figdir**
     - str
     - Set path for storing figures
