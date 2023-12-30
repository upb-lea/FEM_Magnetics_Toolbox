User Guide: How to create a FEMMT model
=======================================

This guide explains how a model can be created in femmt and provides all
the necessary information to work with femmt. Many examples for femmt
models can be found in the example folder.

Working directory
--------------------

Every femmt model has a working directory which can be set when creating
an instance of the femmt base class called ``MagneticComponent``. When
running the simulation many files will be created in the working
directory including the model, mesh and multiple result files. It also
contains the electro_magnetic_log.json which the most important
simulation results (e.g. losses, inductance, ...).

Besides the working directory a ``MagneticComponent`` also needs a
``ComponentType``. Currently this can be '``Inductor``', '``Transformer``' or
'``IntegratedTransformer``'.

This results to the following code:

.. code:: python

   import femmt as fmt

   geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.FreqDomain,
                                component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

The ``simulation_type`` specifies the type of simulation to be performed.

- If set to ``FreqDomain``, indicating a frequency domain simulation.
- If set to ``TimeDomain``, indicating a time domain simulations.

The ``Verbosity`` controls the level of detail in the output.

- If set to ``ToConsole``, all outputs are printed in the terminal.
- If set to ``ToFile``, all outputs are written to files.
- If set to ``Silent``, no command line outputs are shown.

This simple feature significantly speeds up simulation time, especially for many automated simulations.

The ``is_gui`` bool asks at first startup for onelab-path. It is a distinction between GUI and command line and is defaults to False in command-line-mode.

Creating a core
------------------

In general, only 2D rotationally symmetric geometries are represented in
FEMMT. Other core shapes must first be converted to a 2D rotationally
symmetric shape. The corresponding values for this (diameter core,
dimensions of the winding window) are taken from the data sheet.
Afterwards, a corresponding geometry is generated automatically.

The following graphics always show only a sectional view of the core
geometry.

.. image:: ../images/geometry_translated.png
	:width: 800

After creating a MagneticComponent a core needs to be created. The core
needs spatial parameters as well as material parameters. The neccessary
spatial parameters are shown in the image below.

.. image:: ../images/geometry_core.png
	:width: 800

Core spatial parameters can be entered manually but FEMMT provides a
database of different practical cores. This database can be accessed
using:

.. code:: python

   core_db = fmt.core-database()["PQ 40/40"]

Now the core object can be created and added to the model (geo object)

.. code:: python

   core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

   core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False,
                    material=mdb.Material.N49, temperature=45, frequency=inductor_frequency,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)
   geo.set_core(core)

Material database
~~~~~~~~~~~~~~~~~

The `material database <https://github.com/upb-lea/materialdatabase>`_ was already introduced in the upper code example with the `material=` parameter. The temperature as well as the frequency are necessary to pick the corresponding data from the datasheet.


Adding air gaps to the core
------------------------------

In the next steps air gaps can be added. Currently it is only possible
to add air gaps in the center leg, there for the 'AirGapLegPosition' is
always 'CenterLeg'. To set the vertical position for a air gap multiple
methods are available:

-  **Center**: The air gap will always be positioned in the center
-  **Percent**: A value between 0 and 100 can be given. Where 0
   represents the bottom end and 100 the top end of the winding window.
-  **Manually**: The specific y coordinate needs to be entered
   manually.

.. image:: ../images/geometry_air_gap.png
	:width: 800

Have a look at the following example on how to create an air gap object
and add it to the model:

.. code:: python

   air_gaps = fmt.AirGaps(method=fmt.AirGapMethod.Percent, core=core)
   air_gaps.add_air_gap(leg_position=fmt.AirGapLegPosition.CenterLeg, height=0.0005, position_value=50)
   geo.set_air_gaps(air_gaps)

Adding an air_gap object is not necessary. If no air gap is needed,
don't add the air gap object to the model.

Set insulation distances
---------------------------

There are multiple insulations implemented in femmt. Some of them are
created as rectangles in the model, some are just adding an offset to
the windings.

Core insulations are the insulations which are created as rectangles in
the model. 4 core insulations will be added: top, left, bottom, right.
The distance of those values can be set with the 'add_core_insulations'
function.

Furthermore there are offset insulations between each turn in the same
winding, a distance between 2 windings in one virtual winding window and
a distance between each virtual winding window. The first two are set
using the 'add_winding_insulations' functions, the last one when
creating such a virtual winding window (vww).

The function '``add_winding_insulations``' therefore contains lists which represent the insulations between turns of the same winding. This does not correspond to
the order conductors are added to the winding! Instead, the winding number is important. The conductors are sorted by ascending winding number.
The lowest winding number therefore is combined with index 0. The second lowest with index 1 and so on.

.. image:: ../images/geometry_insulation.png
	:width: 800

This is how to create an insulation object and add certain insulations:

.. code:: python

   insulation = fmt.Insulation(flag_insulation=True)
   insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
   insulation.add_winding_insulations([[0.0002, 0.001],[0.001, 0.0002]])
   geo.set_insulation(insulation)

The spatial parameters for the insulation, as well as for every other
function in FEMMT, are always in SI-Units, in this case metres.

Add windings to the winding window
-------------------------------------

In order to understand the way winding windows work in femmt, the
concept of virtual winding windows must be explained:

Virtual Winding Windows
~~~~~~~~~~~~~~~~~~~~~~~

For every femmt model there is always one winding window, which is a 2D
representation of the 3D rotated winding window. This winding window can
be split into multiple virtual winding windows which are used to draw
the conductors. There are multiple ways to split a winding window:

-  **NoSplit**: Only 1 virtual winding window will be returned and it
   has the same size as the real winding window.
-  **HorizontalSplit**: 2 virtual winding windows will be returned, one
   for the top and one for the bottom part. The height of the splitting
   line can be set using a horizontal_split_factor (value between 0 and
   1)
-  **VerticalSplit**: 2 virtual winding windows will be returned, one
   for the left and one for the right part. The radius (x-coordinate) of
   the splitting line can be set using a vertical_split_factor (value
   between 0 and 1)
-  **HorizontalAndVerticalSplit**: 4 virtual winding windows are
   returned. One for each corner (in the following order): top_left,
   top_right, bottom_left, bottom_right. In this case the horizontal and
   vertical split factors can be used to set the sizes of each grid
   cell.

.. image:: ../images/geometry_winding_windows.png
	:width: 800
	
In addition to that 2 virtual winding windows can be combined to one
(this is not possible for (top_left, bottom_right) or (top_right,
bottom_left) combinations). This is done using the combine_vww()
function of the WindingWindow class.

Each virtual winding window can be filled with either one single winding
or one interleaved winding.

A winding window with only one virtual winding window can be craeted
like this:

.. code:: python

   winding_window = fmt.WindingWindow(core, insulation)
   vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

Winding types and winding schemes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table gives an overview of the different winding types,
winding schemes and conductor arrangements:

+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
| **WindingType**  | **ConductorType**      | **WindingScheme**         | **ConductorArrangement**  | **WrapParaType**  | **status**       | **description**                                                                    |
+==================+========================+===========================+===========================+===================+==================+====================================================================================+
| Interleaved      |                        |                           |                           |                   |                  | Always needs 2 conductors                                                          |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  | RoundSolid, RoundLitz  |                           |                           |                   |                  |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | Bifilar                   |                           |                   | not implemented  | TODO                                                                               |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | VerticalAlternating       |                           |                   | not implemented  | primary and secondary winding are interleaved vertically (rows)                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | HorizontalAlternating     |                           |                   | implemented      | primary and secondary winding are interleaved horizontally (cols)                  |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | VerticalStacked           |                           |                   | implemented      | primary winding is drawn bottom to top, seoncdary winmding is drawn top to bottom  |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           | Square                    |                   | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           | Hexagonal                 |                   | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  | RectangularSolid       |                           |                           |                   | not implemented  |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
| Single           |                        |                           |                           |                   |                  | Always needs 1 conductor                                                           |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  | RoundSolid, RoundLitz  |                           |                           |                   |                  |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | None                      |                           |                   | implemented      |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           | Square                    |                   | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           | Square full width         |                   | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           | Hexagonal                 |                   | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  | RectangularSolid       |                           |                           |                   |                  |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | Full                      |                           |                   | implemented      | whole virtual winding window contains is filled with one turn                      |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | FoilHorizontal (stacked)  |                           |                   | implemented      | foils are very long (x-axis) and drawn along y-axis                                |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | Square full width         |                           |                   | not implemented  | foils are drawn along x-axis first and then along y-axis                           |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        | FoilVertical              |                           |                   | implemented      | foils are very tall (y-axis) and drawn along x-axis                                |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           |                           | Fixed Thickness   | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+
|                  |                        |                           |                           | Interpolate       | "                |                                                                                    |
+------------------+------------------------+---------------------------+---------------------------+-------------------+------------------+------------------------------------------------------------------------------------+

ConductorArrangement
^^^^^^^^^^^^^^^^^^^^

-  **Square**: conductors are set in next to each other in a grid
-  **Hexagonal**: similar to square but in this case the conductors
   frpmo the next column slips in the free space between two conductors
   from the first column
-  **Square full width**: conducors are first drawn along x-axis and
   then y-axis

WrapParaType
^^^^^^^^^^^^

-  **Fixed thickness**: TODO
-  **Interpolate**: TODO

Images for the possible winding types can be found :ref:`here <winding_types>`.


Add conductors
-----------------

When creating an instance of the class Conductor a winding number and a
conductivity needs to be given:

The winding number represents the index of the winding (e.g. primary->1,
secondary->2, tertiary->3). As an example: When starting a simulation on
a transformer a current needs to be given, this is done in a list. The
first index of the current's list will be set to the winding with the
lowest winding number, the second index of the list to the winding with
the second lowest winding number and so on.

The conductivity can be set using the Conductivity enum where one of two
possible materials need to be selected:

-  **Copper**
-  **Aluminium**

After creating an conductor object it is necessary to add a conductor to
it. As already shown in the winding types table 3 different conducors
can be set:

-  **RoundSolid**
-  **RoundLitz**
-  **RectangularSolid**

To create a conductor have a look at the following code example:

.. code:: python

   winding1 = fmt.Conductor(winding_number=0, conductivity=fmt.Conductivity.Copper)
   winding1.set_solid_round_conductor(conductor_radius=0.0011, conductor_arrangement=fmt.ConductorArrangement.Square)

Add conductors to virtual winding windows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now the conductors need to be added to the virtual winding windows with
the corresponding winding type and winding scheme. In this case the
set_winding() or set_interleaved_winding() function needs to be called.
In the set_interleaved_winding() function an insulation distance can
also be set. This value represents the distance between conductors from
the primary and secondary side.

.. code:: python

   vww.set_winding(conductor=winding1, turns=9, winding_scheme=None)

If you have a look at the winding types and windng schemes table a
winding scheme is not needed when creating a round solid conductor in
single winding. Therefore the value is set to None.

Now before simulating the winding window needs to be added to the model
as well:

.. code:: python

   geo.set_winding_window(winding_window=winding_window)

Create model and start simulation
------------------------------------

After every needed component is added to the model the model can be
created. This is done using the create_model() function. The frequency
is needed there because of the mesh which is adapted according to the
skin depth. In addition to that a boolean can be given to show the model
after creation (in gmsh).

The last step is to run a simulation using single_simulation(), which
needs the frequency, currents (and phase if transformer is set) as
parameters.

.. code:: python

   geo.create_model(freq=100000, visualize_before=True, save_png=False)
   geo.single_simulation(freq=100000, current=[4.5], show_results=True)

The results should look like this:

[Optional] Create thermal simulation
---------------------------------------

After running the electromagnetic simulation it is possible to use the
simulation results and the created model and start a thermal simulation.
The thermal simulation will add a case surrounding the previous created
model. At the edge of this case the boundary condition is applied and
the thermal conductivity as well as the dimensions of the case can be
choosen freely. This case is split into 5 parts: top, top right, right,
bot right, bot. For each region a different thermal conductivity and
boundary condition can be set. In order to run thermal a thermal
simulation different values are needed:

-  thermal conductivity dict: A dictionary containing thermal
   conductivities for each region. The regions are: air, core, winding,
   air_gaps, insulation, case (which is split in top, top_right, right,
   bot_right, bot
-  case gap values: Set the size of the surrounding case
-  boundary temperatures dict: The temperatures which will be applied at
   the edge of the case (dirichlet boundary condition)
-  boundary flags: By disabling a specific boundary its condition can be
   set to a neumann boundary condition ignoring the temperature
   parameter
   
.. image:: ../images/geometry_thermal.png
	:width: 800

Have a look at this example on how to set the parameters since the
dictionary keywords are important to write correctly:

.. code:: python

   thermal_conductivity_dict = {
           "air": 0.0263,
           "case": {
               "top": 0.122,
               "top_right": 0.122,
               "right": 0.122,
               "bot_right": 0.122,
               "bot": 0.122
           },
           "core": 5,
           "winding": 400,
           "air_gaps": 180,
           "insulation": 0.42
   }

   case_gap_top = 0.002
   case_gap_right = 0.0025
   case_gap_bot = 0.002

   boundary_temperatures = {
       "value_boundary_top": 20,
       "value_boundary_top_right": 20,
       "value_boundary_right_top": 20,
       "value_boundary_right": 20,
       "value_boundary_right_bottom": 20,
       "value_boundary_bottom_right": 20,
       "value_boundary_bottom": 20
   }

   boundary_flags = {
       "flag_boundary_top": 0,
       "flag_boundary_top_right": 0,
       "flag_boundary_right_top": 1,
       "flag_boundary_right": 1,
       "flag_boundary_right_bottom": 1,
       "flag_boundary_bottom_right": 1,
       "flag_boundary_bottom": 1
   }

In the boundary_flags dictionary 2 flags are set to 0 which means there
will be a neumann boundary applied. Please have a look at the picture
above which shows the current selected boundaries.

In the following table a possible set of thermal conductivities can be
found: 

+-----------------------------+--------------------------+
| **Material**                | **Thermal conductivity** |
+-----------------------------+--------------------------+
| air (background)            | 0.0263                   |
+-----------------------------+--------------------------+
| epoxy resign (used in case) | 1.54                     |
+-----------------------------+--------------------------+
| ferrite (core)              | 5                        |
+-----------------------------+--------------------------+
| copper (winding)            | 400                      |
+-----------------------------+--------------------------+
| aluminiumnitride (air gaps) | 180                      |
+-----------------------------+--------------------------+
| polyethylen (insulation)    | 0.42                     |
+-----------------------------+--------------------------+

The thermal simulation will solve the stationary heat equation and since
no convection is considered every material is assumed to be solid. Now
the simulation can be run:

.. code:: python

   geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot, True, True)

The following image shows the simulation results:

.. image:: ../images/user_guide_example_thermal_simulation.png
	:width: 350


