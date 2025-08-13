import os
import femmt as fmt
import logging
import numpy as np
import materialdatabase as mdb

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

# Working directory can be set arbitrarily
working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
if not os.path.exists(working_directory):
    os.mkdir(working_directory)

geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Inductor,
                            working_directory=working_directory, onelab_verbosity=fmt.Verbosity.Silent, is_gui=True)

# 2. set core parameters
core_db = fmt.core_database()["PQ 40/40"]
core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                window_w=core_db["window_w"],
                                                window_h=core_db["window_h"],
                                                core_h=core_db["core_h"])
inductor_frequency = 270000
core = fmt.Core(core_type=fmt.CoreType.Single,
                core_dimensions=core_dimensions,
                material=mdb.Material.N49, temperature=45,
                permeability_datasource=fmt.MaterialDataSource.Custom,
                mu_r_abs=3000, phi_mu_deg=0,
                permittivity_datasource=fmt.MaterialDataSource.Custom,
                mdb_verbosity=fmt.Verbosity.Silent,
                sigma=1)
geo.set_core(core)

# 3. set air gap parameters
air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
geo.set_air_gaps(air_gaps)

# 4. set insulations
insulation = fmt.Insulation(flag_insulation=True)
insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
insulation.add_winding_insulations([[0.0005]])
geo.set_insulation(insulation)

# 5. create winding window and virtual winding windows (vww)
winding_window = fmt.WindingWindow(core, insulation)
vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

# 6. create conductor and set parameters: use solid wires
winding = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=45)
winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
winding.parallel = False

# 7. add conductor to vww and add winding window to MagneticComponent
vww.set_winding(winding, 7, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
geo.set_winding_windows([winding_window])

# 8. create the model
geo.create_model(freq=270000, pre_visualize_geometry=False, save_png=False)

# 6.a. start simulation
# time value list
t = np.linspace(0, 1 / inductor_frequency, 5)
t_list = [float(x) for x in t.tolist()]
#  Current values list
current_values = 4.5 * np.cos(2 * np.pi * inductor_frequency * t)
current_values_list = current_values.tolist()

geo.time_domain_simulation(current_period_vec=[current_values_list],
                           time_period_vec=t_list,
                           number_of_periods=1,
                           show_fem_simulation_results=False,
                           show_rolling_average=False, )
