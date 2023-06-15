import femmt as fmt
import os


example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

# Example for a transformer with multiple virtual winding windows.
working_directory = os.path.join(example_results_folder, "n-horizontal-split")
if not os.path.exists(working_directory):
    os.mkdir(working_directory)

# 1. chose simulation type
geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                            silent=False)

# 2. set core parameters
core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=0.018, window_w=0.02, core_inner_diameter=0.01)
core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
                permeability_datasource=fmt.MaterialDataSource.Custom,
                permittivity_datasource=fmt.MaterialDataSource.Custom)
geo.set_core(core)

# 3. set air gap parameters
air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
geo.set_air_gaps(air_gaps)

# 4. set insulation
insulation = fmt.Insulation()
insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
# insulation.add_winding_insulations([0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002], 0.0005)
iso_self = 0.0001
iso_against = 0.0002
insulation.add_winding_insulations(
    [[iso_self, iso_against, iso_against, iso_against, iso_against, iso_against, iso_against],
     [iso_against, iso_self, iso_against, iso_against, iso_against, iso_against, iso_against],
     [iso_against, iso_against, iso_self, iso_against, iso_against, iso_against, iso_against],
     [iso_against, iso_against, iso_against, iso_self, iso_against, iso_against, iso_against],
     [iso_against, iso_against, iso_against, iso_against, iso_self, iso_against, iso_against],
     [iso_against, iso_against, iso_against, iso_against, iso_against, iso_self, iso_against],
     [iso_against, iso_against, iso_against, iso_against, iso_against, iso_against, iso_self]])
geo.set_insulation(insulation)

# 5. create winding window and virtual winding windows (vww)
winding_window = fmt.WindingWindow(core, insulation)
cells = winding_window.NHorizontalSplit([1/7, 2/7, 3/7, 4/7, 5/7, 6/7])
# 6. create conductors and set parameters

winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
winding1.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
winding3.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
winding4.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

winding5 = fmt.Conductor(4, fmt.Conductivity.Copper)
winding5.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

winding6 = fmt.Conductor(5, fmt.Conductivity.Copper)
winding6.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

winding7 = fmt.Conductor(6, fmt.Conductivity.Copper)
winding7.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)


# 5. create winding window and virtual winding windows (vww)
# winding_window = fmt.WindingWindow(core, insulation)
# cells = winding_window.NCellsSplit(0, [1/6, 2/6, 3/6, 4/6, 5/6])
# 7. add conductor to vww and add winding window to MagneticComponent
cells[0].set_winding(winding1, 6, fmt.WindingType.Single)
cells[1].set_winding(winding2, 5, fmt.WindingType.Single)
cells[2].set_winding(winding3, 5, fmt.WindingType.Single)
cells[3].set_winding(winding4, 5, fmt.WindingType.Single)
cells[4].set_winding(winding5, 5, fmt.WindingType.Single)
cells[5].set_winding(winding6, 5, fmt.WindingType.Single)
cells[6].set_winding(winding6, 3, fmt.WindingType.Single)
geo.set_winding_windows([winding_window])
# 8. start simulation with given frequency, currents and phases
geo.create_model(freq=250000, pre_visualize_geometry=True)
# geo.single_simulation(freq=250000, current=[4, 4, 4, 4, 4, 4, 4],
#                       phi_deg=[0, 180, 0, 180, 0, 180, 90], show_fem_simulation_results=True)

# read inductances
geo.get_inductances(I0=1, op_frequency=250000)
# Reference simulation using FEMM
# geo.femm_reference(freq=250000, current=[4, 4, 4, 4, 4, 4, 4, 4, 4, 4], sign=[1, -1, -1, -1, -1, -1, -1, -1, -1, -1], non_visualize=0)
