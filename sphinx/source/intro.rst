.. sectnum::


.. include:: ../../README.rst


FEMMT class and function documentation
=======================================
The ``MagneticComponent`` class
---------------------------------------
.. currentmodule:: femmt.MagneticComponent

.. autoclass:: MagneticComponent
   :members: set_core, set_air_gaps, set_insulation, set_stray_path, set_winding_window, create_model, single_simulation, excitation_sweep, mesh, thermal_simulation, femm_reference, femm_thermal_validation
   :special-members: __init__

The ``Core`` class
--------------------------------------
.. autoclass:: femmt.Model.Core
    :special-members: __init__
	
The ``AirGaps`` class
--------------------------------------
.. autoclass:: femmt.Model.AirGaps
    :members: add_air_gap
    :special-members: __init__
	
The ``Insulation`` class
--------------------------------------
.. autoclass:: femmt.Model.Insulation
    :members: add_winding_insulations, add_core_insulations
    :special-members: __init__
	
The ``Conductor`` class
--------------------------------------
.. autoclass:: femmt.Model.Conductor
    :members: set_rectangular_conductor, set_solid_round_conductor, set_litz_round_conductor
    :special-members: __init__
	
The ``WindingWindow`` class
--------------------------------------
.. autoclass:: femmt.Model.WindingWindow
    :members: split_window, combine_vww
    :special-members: __init__
	
The ``VirtualWindingWindow`` class
--------------------------------------
.. autoclass:: femmt.Model.VirtualWindingWindow
    :members: set_winding, set_interleaved_winding
    :special-members: __init__
	
The ``LogParser`` class
-----------------------------------
.. autoclass:: femmt.FEMMTLogParser
    :members: plot_frequency_sweep_losses, plot_frequency_sweep_winding_params
    :special-members: __init__
	
.. autoclass:: femmt.FileData
    :members:
    :undoc-members:
	
.. autoclass:: femmt.SweepData
    :members:
    :undoc-members:
	
.. autoclass:: femmt.WindingData
    :members:
    :undoc-members:

``Enumerations``
---------------------------------
.. automodule:: femmt.Enumerations
    :members:
    :undoc-members:

Helper functions
---------------------------------
.. automodule:: femmt.Functions
    :members: core_database, litz_database, wire_material_database

