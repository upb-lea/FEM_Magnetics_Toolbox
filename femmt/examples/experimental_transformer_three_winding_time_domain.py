"""Basic example to show how to simulate a three winding transformer in time domain."""
import femmt as fmt
import materialdatabase as mdb
import os
# from matplotlib import pyplot as plt
import numpy as np
import logging

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def basic_example_transformer_three_windings_time_domain(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Demonstrate how to simulate a three winding transformer in time domain.

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

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain,
                                component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                is_gui=is_test)

    # This line is for automated pytest running on github only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=0.06, window_w=0.03, core_inner_diameter=0.015, core_h=0.08)
    # core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
    #                 permeability_datasource=fmt.MaterialDataSource.Custom,
    #                 permittivity_datasource=fmt.MaterialDataSource.Custom)
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

    # 4. set insulation
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([[0.0002, 0.0004, 0.0004],
                                        [0.0004, 0.0002, 0.0004],
                                        [0.0004, 0.0004, 0.0002]])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    top_left, top_right, bot_left, bot_right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalAndVerticalSplit, horizontal_split_factor=0.4)
    top_left = winding_window.combine_vww(top_left, bot_left)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
    winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)
    winding2 = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
    winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)
    winding3 = fmt.Conductor(2, fmt.ConductorMaterial.Copper)
    winding3.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    top_left.set_winding(winding1, 7, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)
    top_right.set_winding(winding2, 8, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)
    bot_right.set_winding(winding3, 9, fmt.WindingType.Single, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)

    geo.set_winding_windows([winding_window])

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, pre_visualize_geometry=show_visual_outputs)

    # t = np.linspace(0, 2, 30) * 1/inductor_frequency
    t = np.linspace(0, 2 / 250000, 5)
    t_list = [float(x) for x in t.tolist()]
    # # Current values
    current_values_1 = 4 * np.cos(2 * np.pi * 250000 * t)
    current_values_2 = 4 * np.cos(2 * np.pi * 250000 * t + np.pi)
    current_values_3 = 4 * np.cos(2 * np.pi * 250000 * t)
    current_values_list_1 = current_values_1.tolist()
    current_values_list_2 = current_values_2.tolist()
    current_values_list_3 = current_values_3.tolist()

    print(t_list)
    print(current_values_list_1)

    # time_list = [0, 2, 4, 6, 8]
    # plt.plot(t_list, current_values_list_1)
    # plt.plot(t_list, current_values_list_2)
    # plt.plot(t_list, current_values_list_3)
    # plt.xlabel('Time (s)')
    # plt.ylabel('Current (A)')
    # plt.title(f'Cos wave: {250000} Hz, {4.5} A amplitude')
    # plt.grid(True)
    # if show_visual_outputs and not is_test:
    #     plt.show()
    # Electromagnetic time-domain simulation
    # The 'current_period_vec' parameter accepts a list of lists, where each sublist represents the current values for a particular winding.
    # The 'time_period_vec' parameter accepts a single list representing the time steps for the simulation; this is common for all windings.
    # The 'number_of_periods' parameter defines the number of periods or the simulation duration time internally.
    # The 'show_rolling_average' parameter is a boolean that determines whether to plot the rolling average of the simulation results.
    # The 'rolling_avg_window_size' parameter is an integer that specifies the window size for calculating the rolling average.
    # It defines the number of data points used in each calculation of the average
    # a too-small window size may not effectively smooth out short-term fluctuations.
    geo.time_domain_simulation(current_period_vec=[current_values_list_1, current_values_list_2, current_values_list_3],
                               time_period_vec=t_list,
                               number_of_periods=1,
                               plot_interpolation=False,
                               show_fem_simulation_results=show_visual_outputs,
                               show_rolling_average=False,
                               rolling_avg_window_size=5)


if __name__ == "__main__":
    basic_example_transformer_three_windings_time_domain(show_visual_outputs=True)
