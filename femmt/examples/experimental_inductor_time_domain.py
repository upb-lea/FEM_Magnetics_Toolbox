"""Basic example to show how to simulate an inductor in time domain."""
import numpy as np
import femmt as fmt
import os
import logging

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def basic_example_inductor_time_domain(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """
    Demonstrate how to simulate an inductor in time domain.

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

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain,
                                component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 270000

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core_material = fmt.LinearComplexCoreMaterial(mu_r_abs=3100,
                                                  phi_mu_deg=0,
                                                  dc_conductivity=1,
                                                  eps_r_abs=0,
                                                  phi_eps_deg=0,
                                                  mdb_verbosity=fmt.Verbosity.Silent)

    core = fmt.Core(material=core_material,
                    core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False)

    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0002, 90)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([[0.0005]], per_layer_of_turns=False)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=45)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
    # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 7, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward, zigzag=False)
    geo.set_winding_windows([winding_window])

    # 8. create the model
    geo.create_model(freq=inductor_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)

    # 6.a. start simulation
    # time value
    t = np.linspace(0, 1 / inductor_frequency, 5)
    t_list = [float(x) for x in t.tolist()]
    # Current values
    current_values = 4.5 * np.cos(2 * np.pi * inductor_frequency * t)
    current_values_list = current_values.tolist()  # Convert numpy array to list
    print(t_list)
    print(current_values_list)
    # plot to see the current
    # plt.plot(t_list, current_values_list)
    # plt.xlabel('Time (s)')
    # plt.ylabel('Current (A)')
    # plt.title(f'Cos wave: {inductor_frequency} Hz, {4.5} A amplitude')
    # plt.grid(True)
    # if show_visual_outputs and not is_test:
    #     plt.show()

    # Electromagnetic time-domain simulation
    #  The 'current_period_vec' parameter accepts a list of lists, where each sublist represents the current values for a particular winding.
    #  The 'time_period_vec' parameter accepts a single list representing the time steps for the simulation; this is common for all windings.
    #  The 'number_of_periods' parameter defines the number of periods or the simulation duration time internally.
    # The 'show_rolling_average' parameter is a boolean that determines whether to plot the rolling average of the simulation results.
    # The 'rolling_avg_window_size' parameter is an integer that specifies the window size for calculating the rolling average.
    # It defines the number of data points used in each calculation of the average
    # a too-small window size may not effectively smooth out short-term fluctuations.
    geo.time_domain_simulation(current_period_vec=[current_values_list],
                               time_period_vec=t_list,
                               number_of_periods=1,
                               plot_interpolation=False,
                               show_fem_simulation_results=show_visual_outputs,
                               show_rolling_average=False,
                               rolling_avg_window_size=5)


if __name__ == "__main__":
    basic_example_inductor_time_domain(show_visual_outputs=True)
