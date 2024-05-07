.. sectnum::

.. include:: readme.rst

.. include:: introduction.rst

.. _user_guide_model_creation:

.. include:: user_guide_model_creation.rst

.. _winding_types:

.. include:: winding_types.rst

FEMMT class and function documentation
=======================================
The ``MagneticComponent`` class
---------------------------------------
.. currentmodule:: femmt.MagneticComponent

.. autoclass:: femmt.MagneticComponent
   :members: set_core, set_air_gaps, set_insulation, set_stray_path, create_model, single_simulation, excitation_sweep, mesh, thermal_simulation, femm_reference, femm_thermal_validation, set_winding_windows
   :special-members: __init__

The ``Core`` class
--------------------------------------
.. autoclass:: femmt.Core
    :special-members: __init__
	
The ``AirGaps`` class
--------------------------------------
.. autoclass:: femmt.AirGaps
    :members: add_air_gap
    :special-members: __init__
	
The ``Insulation`` class
--------------------------------------
.. autoclass:: femmt.Insulation
    :members: add_winding_insulations, add_core_insulations
    :special-members: __init__, add_core_insulations, add_winding_insulations
	
The ``Conductor`` class
--------------------------------------
.. autoclass:: femmt.Conductor
    :members: set_rectangular_conductor, set_solid_round_conductor, set_litz_round_conductor
    :special-members: __init__
	
The ``WindingWindow`` class
--------------------------------------
.. autoclass:: femmt.WindingWindow
    :members: split_window, combine_vww
    :special-members: __init__, split_window
	
The ``VirtualWindingWindow`` class
--------------------------------------
.. autoclass:: femmt.VirtualWindingWindow
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
.. currentmodule:: femmt
.. autoclass:: femmt.enumerations

Helper functions
---------------------------------
.. automodule:: femmt.functions
    :members:

Model helper functions
---------------------------------
.. automodule:: femmt.functions_model
    :members:

Reluctance model helper functions
---------------------------------
.. automodule:: femmt.reluctance
    :members:

Topology helper functions
---------------------------------
.. automodule:: femmt.functions_topologies
    :members:

Constants
---------------------------------
.. automodule:: femmt.constants
    :members:

.. include:: developer_notes.rst