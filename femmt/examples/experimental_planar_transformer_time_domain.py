"""
Basic example to show how to simulate a PCB transformer with interleaved foil winding.

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
import numpy as np
import matplotlib.pyplot as plt

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


def basic_experimental_planar_transformer_interleaved(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Run the example code for the PCB transformer.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                is_gui=is_test, simulation_type=fmt.SimulationType.TimeDomain, visualization_mode=fmt.VisualizationMode.Stream)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=5e-3,
                                                    window_w=6e-3,
                                                    window_h=2e-3,
                                                    core_h=8e-3)
    core_material = fmt.LinearComplexCoreMaterial(mu_r_abs=3100,
                                                  phi_mu_deg=12,
                                                  dc_conductivity=0.6,
                                                  eps_r_abs=0,
                                                  phi_eps_deg=0)
    core = fmt.Core(material=core_material,
                    core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.0005, 0.0005, 0.0005, 0.0005)
    insulation.add_winding_insulations([[1.15e-4, 1.15e-4],
                                        [1.15e-4, 1.15e-4]])
    insulation.add_insulation_between_layers(thickness=0.01e-3)
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    # vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    # vww_bot, vww_top = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, split_distance=0.0001)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
    winding1.set_rectangular_conductor(thickness=35e-6, width=1.524e-3)

    winding2 = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
    winding2.set_rectangular_conductor(thickness=35e-6, width=4.8e-3)
    winding2.parallel = True

    vww.set_interleaved_winding(winding1, 5, winding2, 2, fmt.InterleavedWindingScheme.VerticalAlternating,
                                foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalUpward, group_size=1)
    geo.set_winding_windows([winding_window])

    # Magnetic Simulation
    # t = np.linspace(0, 2, 30) * 1/inductor_frequency
    t = np.linspace(0, 2 / 1e3, 5)
    t_list = [float(x) for x in t.tolist()]
    # Current values
    current_values_1 = 2 * np.cos(2 * np.pi * 1e3 * t)
    current_values_2 = 2 * np.cos(2 * np.pi * 1e3 * t + np.pi)
    current_values_list_1 = current_values_1.tolist()
    current_values_list_2 = current_values_2.tolist()

    print(len(t_list))
    print(len(current_values_list_1))
    print(len(current_values_list_2))
    time_list = [0, 2, 4, 6, 8]
    plt.plot(t_list, current_values_list_1)
    plt.plot(t_list, current_values_list_2)
    plt.xlabel('Time (s)')
    plt.ylabel('Current (A)')
    plt.title(f'Cos wave: {1e3} Hz, {2} A amplitude')
    plt.grid(True)
    if show_visual_outputs and not is_test:
        plt.show()
    # Create Model
    geo.create_model(freq=1e3, pre_visualize_geometry=show_visual_outputs, save_png=False)


    geo.time_domain_simulation(current_period_vec=[current_values_list_1, current_values_list_2],
                               time_period_vec=t_list,
                               number_of_periods=1,
                               plot_interpolation=False,
                               show_fem_simulation_results=show_visual_outputs,
                               show_rolling_average=False,
                               rolling_avg_window_size=5)


if __name__ == "__main__":
    basic_experimental_planar_transformer_interleaved(show_visual_outputs=True)
