import femmt as fmt
import os

def basic_example_transformer(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):


    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)



    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, "transformer")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                silent=False, is_gui=is_test)

    # This line is for automated pytest running on github only. Please ignore this line!
    if onelab_folder is not None: geo.file_data.onelab_folder_path = onelab_folder

    # 2. set core parameters
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015, window_w=0.012, window_h=0.0295)
    # core = fmt.Core(core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12, sigma=1.2,
    #                 permeability_datasource=fmt.MaterialDataSource.Custom,
    #                 permittivity_datasource=fmt.MaterialDataSource.Custom)
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material="N49", temperature=45, frequency=200000,
                    # permeability_datasource="manufacturer_datasheet",
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    # permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    mu_r_abs=3000, phi_mu_deg=0,
                    permittivity_datasource=fmt.MaterialDataSource.Custom,
                    # permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    sigma=1)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([[0.0002, 0.001],
                                        [0.001, 0.0002]])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    bot, top = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, split_distance=0.001)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    # winding1.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    # winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    # winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)
    winding2.parallel = False
    # winding2.set_litz_round_conductor(0.0011, 50, 0.00011, None, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    bot.set_winding(winding1, 10, None)
    top.set_winding(winding2, 12, None)
    geo.set_winding_windows([winding_window])

    from matplotlib import pyplot as plt
    import numpy as np
    # t = np.linspace(0, 2, 30) * 1/inductor_frequency
    t = np.linspace(0, 2 / 200000, 100)
    t_list = [float(x) for x in t.tolist()]
    # # Current values
    current_values_1 = 2 * np.cos(2 * np.pi * 200000 * t)
    current_values_2 = 2 * np.cos(2* np.pi * 200000 * t + np.pi)
    current_values_list_1 = current_values_1.tolist()
    current_values_list_2 = current_values_2.tolist()

    print(t_list)
    print(current_values_list_1)


    # time_list = [0, 2, 4, 6, 8]
    plt.plot(t_list, current_values_list_1)
    plt.plot(t_list, current_values_list_2)
    plt.xlabel('Time (s)')
    plt.ylabel('Current (A)')
    plt.title(f'Cos wave: {200000} Hz, {4.5} A amplitude')
    plt.grid(True)
    plt.show()

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=200000, pre_visualize_geometry=show_visual_outputs)
    # Electromagnetic time-domain simulation
    #  The 'current' parameter accepts a list of lists, where each sublist represents the current values for a particular winding.
    #  The 'time' parameter accepts a single list representing the time steps for the simulation; this is common for all windings.
    #  The 'time_period' should be always defind as 1/f.. It's introduced to distinguish it from 'timemax'.
    #  The 'initial_time' parameter defines the starting point of the simulation in seconds.
    #  The 'timemax' parameter defines the end point of the simulation in seconds. It can be set to a value like 2*time_period, depending on the user's needs.
    #  The 'NbSteps' parameter represents the total number of time steps within the provided time period for the simulation.
    #  The 'delta_time' is the interval between each time step in seconds, often defined as T/NbSteps.
    # Smaller time steps can capture rapid changes in the system with higher precision, but they also increase the computational load.
    # Note: If 'timemax' is defined as 2 * time_period, it’s essential to adjust 'delta_time' or 'NbSteps' accordingly.
    # For instance, 'delta_time' can be adapted as (2 * time_period / NbSteps) to maintain the simulation's resolution.
    # Failing to do so, the solver will automatically double 'NbSteps', potentially increasing the computational load.
    geo.time_domain_simulation(freq=200000,
                               current=[current_values_list_1, current_values_list_2],
                               time=t_list,
                               time_period=1 / 200000,
                               initial_time=0,
                               timemax=1 / 200000,
                               NbSteps=100,
                               delta_time=((1 / 200000) / 100),
                               plot_interpolation=False,
                               show_fem_simulation_results=True)




if __name__ == "__main__":
    basic_example_transformer(show_visual_outputs=True)