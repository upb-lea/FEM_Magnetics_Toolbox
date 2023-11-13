import numpy as np

import femmt as fmt
import materialdatabase as mdb
import os

def basic_example_inductor(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):



    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, "inductor")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

    # This line is for automated pytest running on github only. Please ignore this line!
    if onelab_folder is not None: geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 270000

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    # core = fmt.Core(core_type=fmt.CoreType.Single,
    #                 core_dimensions=core_dimensions,
    #                 detailed_core_model=False,
    #                 material=mdb.Material.N49, temperature=45, frequency=inductor_frequency,
    #                 # permeability_datasource="manufacturer_datasheet",
    #                 permeability_datasource=fmt.MaterialDataSource.Measurement,
    #                 permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
    #                 permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
    #                 permittivity_datasource=fmt.MaterialDataSource.Measurement,
    #                 permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
    #                 permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)

    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material=mdb.Material.N49, temperature=45, frequency=inductor_frequency,
                    permeability_datasource=fmt.MaterialDataSource.Custom,
                    mu_r_abs=3000, phi_mu_deg=0,
                    permittivity_datasource=fmt.MaterialDataSource.Custom,
                    mdb_verbosity=fmt.Verbosity.Silent,
                    sigma=1)
    # mu_rel=3000, phi_mu_deg=10,
    # sigma=0.5)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0002, 90)
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
    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=45)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
    # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
    # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 7, None)
    geo.set_winding_windows([winding_window])

    # 8. create the model
    geo.create_model(freq=inductor_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)

    # 6.a. start simulation
    from matplotlib import pyplot as plt
    # time value
    t = np.linspace(0, 2 / inductor_frequency, 50)
    # t_list = t.tolist()
    t_list = [float(x) for x in t.tolist()]
    # # Current values
    current_values = 4.5 * np.cos(2 * np.pi * inductor_frequency * t)
    current_values_list = current_values.tolist()  # Convert numpy array to list
    print(t_list)
    print(current_values_list)
    # plot to see the current
    plt.plot(t_list, current_values_list)
    plt.xlabel('Time (s)')
    plt.ylabel('Current (A)')
    plt.title(f'Cos wave: {inductor_frequency} Hz, {4.5} A amplitude')
    plt.grid(True)
    plt.show()

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
    # The 'show_rolling_average' parameter is a boolean that determines whether to plot the rolling average of the simulation results.
    #
    # The 'rolling_avg_window_size' parameter is an integer that specifies the window size for calculating the rolling average.
    # It defines the number of data points used in each calculation of the average
    # a too-small window size may not effectively smooth out short-term fluctuations.
    #
    # The rolling average is particularly useful for analyzing data that contains noise or fluctuations, allowing for a clearer
    # visualization of the underlying trends or patterns in the data. When 'show_rolling_average' is set to True, the rolling
    # average is computed for each result file, and plots displaying both the original data and the rolling average are shown
    geo.time_domain_simulation(freq=inductor_frequency,
                               current=[current_values_list],
                               time=t_list,
                               time_period=1 / inductor_frequency, #hide it
                               initial_time=0,
                               timemax=2 / inductor_frequency,
                               NbSteps=50,
                               delta_time=(2 / inductor_frequency) / 50, #hide it
                               plot_interpolation=False,
                               show_fem_simulation_results=True,
                               show_rolling_average=True,
                               rolling_avg_window_size=50)

if __name__ == "__main__":
    basic_example_inductor(show_visual_outputs=True)

