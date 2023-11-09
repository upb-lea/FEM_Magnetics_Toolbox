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
geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory, verbosity=fmt.Verbosity.ToConsole)

# 2. set core parameters
core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=16.1e-3, window_w=(22.5-12)/2*1e-3, core_inner_diameter=12e-3, core_h=22e-3)
core = fmt.Core(core_dimensions=core_dimensions, material=fmt.Material.N95, temperature=60, frequency=100000,
                    # permeability_datasource="manufacturer_datasheet",
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK)
geo.set_core(core)

# 3. set air gap parameters
air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.00016, 50)
geo.set_air_gaps(air_gaps)

# 4. set insulation
insulation = fmt.Insulation()
insulation.add_core_insulations(0.0008, 0.0008, 0.001, 0.0001)
# insulation.add_winding_insulations([0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002], 0.0005)
iso_self = 0.0001
iso_against = 0.0002
insulation.add_winding_insulations(
    [[iso_self, iso_against, iso_against, iso_against, iso_against],
     [iso_against, iso_self, iso_against, iso_against, iso_against],
     [iso_against, iso_against, iso_self, iso_against, iso_against],
     [iso_against, iso_against, iso_self, iso_against, iso_against],
     [iso_against, iso_against, iso_against, iso_against, iso_self]])
geo.set_insulation(insulation)

# 5. create winding window and virtual winding windows (vww)
winding_window = fmt.WindingWindow(core, insulation)
cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.48, 0.75], vertical_split_factors=[None, [0.5, 0.85], None])

# 6. create windings and assign conductors
winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
winding1.set_litz_round_conductor(0.85e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
winding2.set_litz_round_conductor(1.0e-3 / 2, 60, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
winding3.set_litz_round_conductor(0.75e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

winding4 = fmt.Conductor(3, fmt.Conductivity.Copper)
winding4.set_litz_round_conductor(0.95e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

winding5 = fmt.Conductor(4, fmt.Conductivity.Copper)
winding5.set_litz_round_conductor(0.75e-3 / 2, 40, 0.1e-3 / 2, None, fmt.ConductorArrangement.Square)

# 7. assign windings to virtual winding windows (cells)
cells[0].set_winding(winding1, 22, fmt.WindingType.Single)
cells[1].set_winding(winding2, 6, fmt.WindingType.Single)
cells[2].set_winding(winding3, 6, fmt.WindingType.Single)
cells[3].set_winding(winding4, 1, fmt.WindingType.Single)
cells[4].set_winding(winding5, 2, fmt.WindingType.Single)
geo.set_winding_windows([winding_window])

# 8. perform an FEM simulation
geo.create_model(freq=100000, pre_visualize_geometry=True)
geo.single_simulation(freq=100000, current=[1.625, 6.9, 4.9, 0.5, 1],
                      phi_deg=[0, 180, 180, 90, 89], show_fem_simulation_results=True)
