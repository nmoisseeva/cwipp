****************
Documentation
****************

Plume Classes
-----------------



.. currentmodule:: cwipp

.. autosummary::
   :nosignatures:
   :recursive:
   :toctree:

   cwipp.Plume
   cwipp.LESplume
   cwipp.MODplume


Plume Base Class Methods
------------------------

Common methods, which apply to all plumes.

.. autosummary::
   :recursive:
   :nosignatures:
   :toctree:

   Plume.get_I
   Plume.get_wf
   Plume.classify


LESplume Class Methods
------------------------

These methods can be applied for plumes with full cross-sectional data available from numerical simulations.

.. autosummary::
   :recursive:
   :nosignatures:
   :toctree:

   LESplume.get_zCL


MODplume Class Methods
------------------------

Methods that apply to plumes without full cross-sectional data available (a.k.a forecast mode).

.. autosummary::
  :recursive:
  :nosignatures:
  :toctree:

  MODplume.iterate
  MODplume.explicit_solution
