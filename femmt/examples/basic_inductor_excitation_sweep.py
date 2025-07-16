"""
Example file how to perform some frequency sweeps on an inductor.

Does not contain a thermal simulation.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os
import logging

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

# Create Object
def basic_example_inductor_excitation_sweep(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Run the example code for the inductor excitation sweep.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    # 0: choose frequencies, amplitude and phases to sweep
    frequencies = [100000, 200000]
    current_amplitudes = [[10], [4]]
    phases = [[0], [179]]

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material="N95", temperature=25, frequency=100000,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup="LEA_LK",
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup="LEA_LK")

    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([[0.0005]], per_layer_of_turns=False)
    insulation.add_turn_insulation([0.25e-5], add_turn_insulations=False)
    insulation.add_insulation_between_layers(add_insulation_material=False, thickness=0.0005)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=100)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 9, None)
    geo.set_winding_windows([winding_window])

    # 8. create the model
    geo.create_model(freq=100000, pre_visualize_geometry=show_visual_outputs, save_png=False)

    # 9. start simulation
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases,
                         show_last_fem_simulation=show_visual_outputs)


if __name__ == "__main__":
    basic_example_inductor_excitation_sweep(show_visual_outputs=True)
