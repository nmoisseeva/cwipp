.. PR_model documentation master file, created by
   sphinx-quickstart on Thu Oct 22 17:36:50 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CWIPP model documentation
====================================

Cross-Wind Integrated Plume Penetration (CWIPP) model is an analytical method for predicting the vertical distribution of smoke in the atmosphere above a wildfire.

This collection of diagnostic routines can be applied as a stand-alone model for real-world plumes, as well as used with numerical output from WRF-SFIRE LES simulations.

Scientific background and details of model development can be found in:
https://acp.copernicus.org/preprints/acp-2020-827/

This documentation page is designed to aid implementation of CWIPP within BlueSky Canada smoke modelling framework with the goal of improving operational forecasts produced for:
https://firesmoke.ca/


.. image:: _static/images/intro.png
   :align: center


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   ./install
   ./usage
   ./reference
   ./LES


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
