****************
Documentation
****************

Plume Classes
-----------------



.. currentmodule:: cwipp

.. autosummary::
   :nosignatures:

   cwipp.Plume
   cwipp.LESplume
   cwipp.MODplume


Plume Base Class Methods
------------------------

Common methods, which apply to all plumes.


.. autosummary::
   :nosignatures:
   :toctree:

   Plume.get_sounding
   Plume.get_wf
   Plume.classify


LESplume Class Methods
------------------------

These methods can be applied for plumes with full cross-sectional data available from numerical simulations.


.. autosummary::
   :nosignatures:
   :toctree:

   LESplume.get_I
   LESplume.get_zCL


MODplume Class Methods
------------------------

Methods that apply to plumes without full cross-sectional data available (a.k.a forecast mode).


.. autosummary::
  :nosignatures:
  :toctree:

  MODplume.iterate
  MODplume.explicit_solution
