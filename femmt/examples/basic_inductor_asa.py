"""
Basic example to show how to simulate an inductor.

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

def basic_example_inductor(onelab_folder: str = None, show_visual_outputs: bool = False, is_test: bool = False):
    """
    Run the example code for the inductor.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    # Added by ASA
    target_inductance = 2e-05
    act_air_gap = 0.0040

    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.FreqDomain, component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 270_000

    # 2.1a set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core_material = fmt.ImportedComplexCoreMaterial(material=fmt.Material.N49,
                                                    temperature=45,
                                                    permeability_datasource=fmt.DataSource.MagNet,
                                                    permittivity_datasource=fmt.DataSource.LEA_MTB)

    # ASA Needs to be calculated while initialisatin of core_material! self.database is a problem  (database is to hide)
    initial_mag_curve = core_material.database.get_initial_magnetization_curve(fmt.Material.N49, 500, 45)

    core = fmt.Core(material=core_material,
                    core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False)

    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, act_air_gap, 50)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0002, 90)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    # it is preferred to assign the exact dimensions of the bobbin for running electrostatic simulations or obtaining the capacitance of the inductor component
    # using the function below
    # bobbin_db = fmt.bobbin_database()["PQ 40/40"]
    # bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
    #                                               bobbin_window_w=bobbin_db["bobbin_window_w"],
    #                                               bobbin_window_h=bobbin_db["bobbin_window_h"],
    #                                               bobbin_h=bobbin_db["bobbin_h"])
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(0.001, 0.001, 0.003, 0.001)
    insulation.add_winding_insulations([[0.0005, 0.0005]], per_layer_of_turns=False)
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
    vww.set_winding(winding, 14, None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.HorizontalRightward_VerticalUpward,
                    zigzag=True)
    geo.set_winding_windows([winding_window])


    inductance_list: list[tuple[float,float,float]] = []

    while True:

        # 8. create the model
        geo.create_model(freq=inductor_frequency, pre_visualize_geometry=show_visual_outputs, save_png=False)
        # 6.a. start simulation
        act_ind, st_ind = geo.single_simulation_with_current_offset(freq=inductor_frequency, current=[4.5], current_offset=4.5, initial_mag_curve=initial_mag_curve,
                                                  plot_interpolation=False, show_fem_simulation_results=show_visual_outputs)

        inductance_list.append ((act_ind, act_air_gap, st_ind))
        ratio = act_ind / target_inductance
        if ratio<0.99 or ratio>1.01:
            if len(inductance_list) > 1:
                # Calculate required gap by L_target = gap_factor/(gap - core_constant)
                # 2 measurments are necessary to calculate gap_factor and core_constant
                delta_l2_l1 = inductance_list[-2][0] - inductance_list[-1][0]
                l1_g1 = inductance_list[-1][0] * inductance_list[-1][1]
                l2_g2 = inductance_list[-2][0] * inductance_list[-2][1]
                core_constant = (l1_g1 - l2_g2) / delta_l2_l1
                gap_factor = inductance_list[-1][0] * (inductance_list[-1][1] + core_constant)
                act_air_gap = gap_factor / target_inductance - core_constant
            else:
                act_air_gap = act_air_gap * ratio
            # Overwrite old airgap object
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, act_air_gap, 50)
            geo.set_air_gaps(air_gaps)
            geo.set_air_gaps(air_gaps)
        else:
            break


    # geo.single_simulation(freq=inductor_frequency, current=[4.5],
    #                       plot_interpolation=False, show_fem_simulation_results=show_visual_outputs)
    # geo.get_inductances(I0=2, op_frequency=20000, skin_mesh_factor=0.5)
    # geo.femm_reference(freq=inductor_frequency, current=[4.5], sign=[1], non_visualize=0)#

    # 6.b. Excitation Sweep Example
    # Perform a sweep using more than one frequency
    # fs = [0, 10000, 30000, 60000, 100000, 150000]
    # amplitude_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
    # phase_list = [[0], [10], [20], [30], [40], [50]]
    # geo.excitation_sweep(frequency_list=fs, current_list_list=amplitude_list, phi_deg_list_list=phase_list)
    # 9. start simulation

    # 7. prepare and start thermal simulation
    # example_thermal_simulation(show_visual_outputs, flag_insulation=True)

    # Extract the capacitance of inductor component
    # geo.get_capacitance_of_inductor_component(show_equivalent_circuit=True)
    # geo.get_inductor_stray_capacitance(show_visual_outputs=True)

def basic_example_inductor_asa(onelab_folder: str = None, show_visual_outputs: bool = False, is_test: bool = False):
    """
    Run the example code for the inductor.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    # Added by ASA
    target_inductance = 2e-05
    act_air_gap = 0.0005


    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.FreqDomain, component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                is_gui=is_test)

    # This line is for automated pytest running on GitHub only. Please ignore this line!
    if onelab_folder is not None:
        geo.file_data.onelab_folder_path = onelab_folder

    inductor_frequency = 270_000

    # 2.1a set core parameters
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core_material = fmt.ImportedComplexCoreMaterial(material=fmt.Material.N49,
                                                    temperature=45,
                                                    permeability_datasource=fmt.DataSource.MagNet,
                                                    permittivity_datasource=fmt.DataSource.LEA_MTB)

    # Calculate Jilles-Atherton initial magnetization curve
    initial_mag_curve = core_material.database.get_initial_magnetization_curve(fmt.Material.N49, 500, 45)

    core = fmt.Core(material=core_material,
                    core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    detailed_core_model=False)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, act_air_gap, 50)

    # 4. set insulation
    # it is preferred to assign the exact dimensions of the bobbin for running electrostatic simulations or obtaining the capacitance of the inductor component
    # using the function below
    # bobbin_db = fmt.bobbin_database()["PQ 40/40"]
    # bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
    #                                               bobbin_window_w=bobbin_db["bobbin_window_w"],
    #                                               bobbin_window_h=bobbin_db["bobbin_window_h"],
    #                                               bobbin_h=bobbin_db["bobbin_h"])
    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(0.001, 0.001, 0.003, 0.001)
    insulation.add_winding_insulations([[0.0005, 0.0005]], per_layer_of_turns=False)

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
    vww.set_winding(winding, 14, None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.HorizontalRightward_VerticalUpward,
                    zigzag=True)

    # Initialize the magnetic component
    geo.initialize_component_parameter("N49", 45, "MagNet", "LEA_MTB", core_dimensions=core_dimensions,
                                       air_gaps=air_gaps, insulation=insulation, winding_windows=[winding_window])

    inductance_list: list[tuple[float,float,float]] = []

    while True:

        # 6.a. start simulation
        act_ind, st_ind = geo.single_simulation_with_current_offset(freq=inductor_frequency, current=[4.5], current_offset=0.2, initial_mag_curve=initial_mag_curve,
                                                  plot_interpolation=False, show_fem_simulation_results=show_visual_outputs)

        # act_ind = float(geo.get_inductance_from_reluctance())
        print(f"static inductance:{st_ind} dynamic inductance:{act_ind}")


        inductance_list.append ((act_ind, act_air_gap, st_ind))
        ratio = act_ind / target_inductance
        if ratio<0.99 or ratio>1.01:
            if len(inductance_list) > 1:
                # Calculate required gap by L_target = gap_factor/(gap - core_constant)
                # 2 measurments are necessary to calculate gap_factor and core_constant
                delta_l2_l1 = inductance_list[-2][0] - inductance_list[-1][0]
                l1_g1 = inductance_list[-1][0] * inductance_list[-1][1]
                l2_g2 = inductance_list[-2][0] * inductance_list[-2][1]
                core_constant = (l1_g1 - l2_g2) / delta_l2_l1
                gap_factor = inductance_list[-1][0] * (inductance_list[-1][1] + core_constant)
                act_air_gap = gap_factor / target_inductance - core_constant
            else:
                act_air_gap = act_air_gap * ratio
            # Overwrite old airgap object
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, act_air_gap, 50)
            geo.set_air_gaps(air_gaps)
            geo.set_air_gaps(air_gaps)
        else:
            break


    # geo.single_simulation(freq=inductor_frequency, current=[4.5],
    #                       plot_interpolation=False, show_fem_simulation_results=show_visual_outputs)
    # geo.get_inductances(I0=2, op_frequency=20000, skin_mesh_factor=0.5)
    # geo.femm_reference(freq=inductor_frequency, current=[4.5], sign=[1], non_visualize=0)#

    # 6.b. Excitation Sweep Example
    # Perform a sweep using more than one frequency
    # fs = [0, 10000, 30000, 60000, 100000, 150000]
    # amplitude_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
    # phase_list = [[0], [10], [20], [30], [40], [50]]
    # geo.excitation_sweep(frequency_list=fs, current_list_list=amplitude_list, phi_deg_list_list=phase_list)
    # 9. start simulation

    # 7. prepare and start thermal simulation
    # example_thermal_simulation(show_visual_outputs, flag_insulation=True)

    # Extract the capacitance of inductor component
    # geo.get_capacitance_of_inductor_component(show_equivalent_circuit=True)
    # geo.get_inductor_stray_capacitance(show_visual_outputs=True)


if __name__ == "__main__":
    # basic_example_inductor(show_visual_outputs=False)
    basic_example_inductor_asa(show_visual_outputs=False)