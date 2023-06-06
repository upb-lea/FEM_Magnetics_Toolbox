import pytest
import os
import json
import femmt as fmt
import deepdiff

def compare_result_logs(first_log_filepath, second_log_filepath, significant_digits=6):
    first_content = None
    second_content = None

    with open(first_log_filepath, "r") as fd:
        first_content = json.loads(fd.read())
        if "date" in first_content["simulation_settings"]:
            del first_content["simulation_settings"]["date"]
        if "working_directory" in first_content["simulation_settings"]:
            del first_content["simulation_settings"]["working_directory"]

    with open(second_log_filepath, "r") as fd:
        second_content = json.loads(fd.read())
        if "date" in second_content["simulation_settings"]:
            del second_content["simulation_settings"]["date"]
        if "working_directory" in second_content["simulation_settings"]:
            del second_content["simulation_settings"]["working_directory"]

    difference = deepdiff.DeepDiff(first_content, second_content, ignore_order=True, significant_digits=significant_digits)
    print(f"{difference = }")


    assert not deepdiff.DeepDiff(first_content, second_content, ignore_order=True, significant_digits=significant_digits)
    # made several tests with the deepdiff command:
    # tried adding not existing keys in one of the dicts: results as expected in an error
    # changed values in very nested dict: results as expected in an error
    # So this command is valid to compare the dicts.

    #return first_content == second_content

def compare_thermal_result_logs(first_log_filepath, second_log_filepath, significant_digits=6):

    with open(first_log_filepath, "r") as fd:
        first_content = json.loads(fd.read())

    with open(second_log_filepath, "r") as fd:
        second_content = json.loads(fd.read())

    difference = deepdiff.DeepDiff(first_content, second_content, ignore_order=True, significant_digits=significant_digits)
    print(f"{difference = }")


    assert not deepdiff.DeepDiff(first_content, second_content, ignore_order=True, significant_digits=significant_digits)
    # made several tests with the deepdiff command:
    # tried adding not existing keys in one of the dicts: results as expected in an error
    # changed values in very nested dict: results as expected in an error
    # So this command is valid to compare the dicts.

    #return first_content == second_content



@pytest.fixture
def temp_folder():
    # Setup temp folder
    temp_folder_path = os.path.join(os.path.dirname(__file__), "temp")

    if not os.path.exists(temp_folder_path):
        os.mkdir(temp_folder_path)

    # Get onelab path
    onelab_path = os.path.join(os.path.dirname(__file__), "..", "..", "onelab")

    # Test
    yield temp_folder_path, onelab_path

@pytest.fixture
def femmt_simulation_inductor_core_material_database(temp_folder):
    temp_folder_path, onelab_folder = temp_folder
    
    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions, material="N95", temperature=25, frequency=100000,
                        permeability_datasource=fmt.MaterialDataSource.ManufacturerDatasheet,
                        permittivity_datasource=fmt.MaterialDataSource.ManufacturerDatasheet)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 9, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=False)

        """
        Currently only the magnetics simulation is tested

        thermal_conductivity_dict = {
                "air": 0.0263,
                "case": 0.3,
                "core": 5,
                "winding": 400,
                "air_gaps": 180
        }
        case_gap_top = 0.0015
        case_gap_right = 0.0025
        case_gap_bot = 0.002
        boundary_temperatures = {
            "value_boundary_top": 293,
            "value_boundary_top_right": 293,
            "value_boundary_right_top": 293,
            "value_boundary_right": 293,
            "value_boundary_right_bottom": 293,
            "value_boundary_bottom_right": 293,
            "value_boundary_bottom": 293
        }
        boundary_flags = {
            "flag_boundary_top": 1,
            "flag_boundary_top_right": 1,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot, show_results=False)
        """
    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")


@pytest.fixture
def femmt_simulation_inductor_core_material_database_measurement(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions, material="N95", temperature=25, frequency=100000,
                        permeability_datasource=fmt.MaterialDataSource.Measurement,
                        permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                        permeability_measurement_setup="LEA_LK",
                        permittivity_datasource=fmt.MaterialDataSource.Measurement,
                        permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                        permittivity_measurement_setup="LEA_LK")
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_solid_round_conductor(conductor_radius=0.0013,
                                          conductor_arrangement=fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 9, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=False)

        """
        Currently only the magnetics simulation is tested

        thermal_conductivity_dict = {
                "air": 0.0263,
                "case": 0.3,
                "core": 5,
                "winding": 400,
                "air_gaps": 180
        }
        case_gap_top = 0.0015
        case_gap_right = 0.0025
        case_gap_bot = 0.002
        boundary_temperatures = {
            "value_boundary_top": 293,
            "value_boundary_top_right": 293,
            "value_boundary_right_top": 293,
            "value_boundary_right": 293,
            "value_boundary_right_bottom": 293,
            "value_boundary_bottom_right": 293,
            "value_boundary_bottom": 293
        }
        boundary_flags = {
            "flag_boundary_top": 1,
            "flag_boundary_top_right": 1,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top, case_gap_right, case_gap_bot, show_results=False)
        """
    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")


@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        mu_r_abs=3000, phi_mu_deg=10, sigma=0.5, permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_solid_round_conductor(conductor_radius=0.0013,
                                          conductor_arrangement=fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 9, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle_litz_wire(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3000, phi_mu_deg=10, sigma=0.5, permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_litz_round_conductor(None, 100, 70e-6, 0.5, fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 9, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Choose wrap para type
        wrap_para_type = fmt.WrapParaType.FixedThickness

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_rectangular_conductor(thickness=1e-3)

        vww.set_winding(winding, 5, fmt.WindingScheme.FoilVertical, wrap_para_type)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[3], show_fem_simulation_results=False)



    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate

    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)
        # Choose wrap para type
        wrap_para_type = fmt.WrapParaType.FixedThickness

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_rectangular_conductor(thickness=1e-3)

        vww.set_winding(winding, 12, fmt.WindingScheme.FoilHorizontal, wrap_para_type)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[3], show_fem_simulation_results=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")


@pytest.fixture
def femmt_simulation_transformer_core_fixed_loss_angle(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder


        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015,
                                                        window_w=0.012,
                                                        window_h=0.0295)



        # 2. set core parameters
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=1.2, permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
        geo.set_air_gaps(air_gaps)

        # 4. set insulation
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0002, 0.0002], [0.0002, 0.0002]])
        geo.set_insulation(insulation)

        # 5. create winding window and virtual winding windows (vww)
        winding_window = fmt.WindingWindow(core, insulation)
        left, right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, 0.0005)

        # 6. create conductors and set parameters
        winding1 = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        winding2 = fmt.Conductor(1, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        # 7. add conductor to vww and add winding window to MagneticComponent
        left.set_winding(winding1, 10, None)
        right.set_winding(winding2, 10, None)
        geo.set_winding_windows([winding_window])

        # 8. start simulation with given frequency, currents and phases
        geo.create_model(freq=250000, pre_visualize_geometry=False)
        geo.single_simulation(freq=250000, current=[4, 4], phi_deg=[0, 178], show_fem_simulation_results=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

@pytest.fixture
def femmt_simulation_transformer_interleaved_core_fixed_loss_angle(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015,
                                                        window_w=0.012,
                                                        window_h=0.0295)



        # 2. set core parameters
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        non_linear=False, sigma=1, re_mu_rel=3200, phi_mu_deg=10, permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
        geo.set_air_gaps(air_gaps)

        # 4. set insulations
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0002, 0.0005], [0.0005, 0.0002]])
        geo.set_insulation(insulation)

        # 5. create winding window and virtual winding windows (vww)
        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        # 6. create conductors and set parameters
        winding1 = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding1.set_solid_round_conductor(0.0011, None)

        winding2 = fmt.Conductor(1, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding2.set_solid_round_conductor(0.0011, None)

        # 7. add conductor to vww and add winding window to MagneticComponent
        vww.set_interleaved_winding(winding1, 21, winding2, 7, fmt.InterleavedWindingScheme.HorizontalAlternating)
        geo.set_winding_windows([winding_window])

        # 8. start simulation with given frequency, currents and phases
        geo.create_model(freq=250000, pre_visualize_geometry=False)
        geo.single_simulation(freq=250000, current=[4, 12], phi_deg=[0, 180], show_fem_simulation_results=False)



    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

@pytest.fixture
def femmt_simulation_transformer_integrated_core_fixed_loss_angle(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                    working_directory=working_directory, silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder


        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.02,
                                                        window_w=0.011,
                                                        window_h=0.03)



        # 2. set core parameters
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom, permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        # 2.1 set stray path parameters
        stray_path = fmt.StrayPath(start_index=0, length=geo.core.core_inner_diameter / 2 + geo.core.window_w - 0.001)
        geo.set_stray_path(stray_path)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 30)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 40)
        geo.set_air_gaps(air_gaps)

        # 4. set insulations
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0002, 0.0005], [0.0005, 0.0002]])
        geo.set_insulation(insulation)

        # 5. create winding window and virtual winding windows (vww)
        # For an integrated transformer it is not necessary to set horizontal and vertical split factors
        # since this is determined by the stray_path
        winding_window = fmt.WindingWindow(core, insulation, stray_path, air_gaps)
        top, bot = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, 0.0001)

        # 6. set conductor parameters
        winding1 = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding1.set_solid_round_conductor(0.0011, None)

        winding2 = fmt.Conductor(1, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding2.set_solid_round_conductor(0.0011, None)

        # 7. add conductor to vww and add winding window to MagneticComponent
        top.set_interleaved_winding(winding1, 3, winding2, 6, fmt.InterleavedWindingScheme.HorizontalAlternating)
        bot.set_interleaved_winding(winding1, 1, winding2, 2, fmt.InterleavedWindingScheme.HorizontalAlternating)
        geo.set_winding_windows([winding_window])

        # 8. start simulation with given frequency, currents and phases
        geo.create_model(freq=250000, pre_visualize_geometry=False)
        geo.single_simulation(freq=250000, current=[8.0, 4.0], phi_deg=[0, 175], show_fem_simulation_results=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    return os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

@pytest.fixture
def thermal_simulation(temp_folder):
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor,
                                    working_directory=working_directory, silent=True, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"])
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=0.,
                        non_linear=False,
                        permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom
                        )
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0001]])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
        winding.set_solid_round_conductor(conductor_radius=0.0015,
                                          conductor_arrangement=fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 8, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)
        geo.single_simulation(freq=100000, current=[3], show_fem_simulation_results=False)

        thermal_conductivity_dict = {
            "air": 1.57,  # potting epoxy resign
            "case": {
                "top": 1.57,
                "top_right": 1.57,
                "right": 1.57,
                "bot_right": 1.57,
                "bot": 1.57
            },
            "core": 5,  # ferrite
            "winding": 400,  # copper
            "air_gaps": 1.57,
            "insulation": 1.57
        }

        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

        boundary_flags = {
            "flag_boundary_top": 1,
            "flag_boundary_top_right": 1,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        # color_scheme = fmt.colors_ba_jonas
        # colors_geometry = fmt.colors_geometry_ba_jonas
        color_scheme = fmt.colors_ba_jonas
        colors_geometry = fmt.colors_geometry_draw_only_lines

        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right,
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False, color_scheme=color_scheme,
                               colors_geometry=colors_geometry)


    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result



def test_inductor_core_material_database(femmt_simulation_inductor_core_material_database):
    """
    The first idea was to compare the simulated meshes with test meshes simulated manually.
    It turns out that the meshes cannot be compared because even slightly differences in the mesh,
    can cause to a test failure, because the meshes are binary files.
    Those differences could even occur when running the simulation on different machines
    -> This was observed when creating a docker image and running the tests.

    Now as an example only the result log will be checked.
    """
    test_result_log = femmt_simulation_inductor_core_material_database

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_inductor_core_material.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_inductor_core_material_database_measurement(femmt_simulation_inductor_core_material_database_measurement):
    """
    The first idea was to compare the simulated meshes with test meshes simulated manually.
    It turns out that the meshes cannot be compared because even slightly differences in the mesh,
    can cause to a test failure, because the meshes are binary files.
    Those differences could even occur when running the simulation on different machines
    -> This was observed when creating a docker image and running the tests.

    Now as an example only the result log will be checked.
    """
    test_result_log = femmt_simulation_inductor_core_material_database_measurement

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_inductor_core_material_measurement.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_inductor_core_fixed_loss_angle(femmt_simulation_inductor_core_fixed_loss_angle):
    test_result_log = femmt_simulation_inductor_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_inductor_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_inductor_core_fixed_loss_angle_litz_wire(femmt_simulation_inductor_core_fixed_loss_angle_litz_wire):
    test_result_log = femmt_simulation_inductor_core_fixed_loss_angle_litz_wire

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_inductor_core_fixed_loss_angle_litz_wire.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)

def test_inductor_core_fixed_loss_angle_foil_vertical(femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical):
    test_result_log = femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_fixed_loss_angle_foil_vertical.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_inductor_core_fixed_loss_angle_foil_horizontal(femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal):
    test_result_log = femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_fixed_loss_angle_foil_horizontal.json")
    compare_result_logs(test_result_log, fixture_result_log)




def test_transformer_core_fixed_loss_angle(femmt_simulation_transformer_core_fixed_loss_angle):
    test_result_log = femmt_simulation_transformer_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_transformer_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_transformer_interleaved_core_fixed_loss_angle(femmt_simulation_transformer_interleaved_core_fixed_loss_angle):
    test_result_log = femmt_simulation_transformer_interleaved_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_transformer_interleaved_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_transformer_integrated_core_fixed_loss_angle(femmt_simulation_transformer_integrated_core_fixed_loss_angle):
    test_result_log = femmt_simulation_transformer_integrated_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_transformer_integrated_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

def test_load_files(temp_folder, femmt_simulation_inductor_core_material_database,
                    femmt_simulation_inductor_core_fixed_loss_angle,
                    femmt_simulation_inductor_core_fixed_loss_angle_litz_wire,
                    femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal,
                    femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical,
                    femmt_simulation_transformer_core_fixed_loss_angle,
                    femmt_simulation_transformer_interleaved_core_fixed_loss_angle,
                    femmt_simulation_transformer_integrated_core_fixed_loss_angle
                    ):
    """
    This function tests if simulations can be set up from a simulation file.
    There is no complete function check, there is just an error-check if the load will fail or not.

    Note: Fixtures are used, to make sure that the latest self-generated log-files can be read.
    It is not okay to use predefined logfiles which where generated once.
    """
    temp_folder_path, onelab_folder = temp_folder

    working_directory = temp_folder_path
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)


        result_log_filepath_list = [femmt_simulation_inductor_core_material_database,
                                    femmt_simulation_inductor_core_fixed_loss_angle,
                                    femmt_simulation_inductor_core_fixed_loss_angle_litz_wire,
                                    femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal,
                                    femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical,
                                    femmt_simulation_transformer_core_fixed_loss_angle,
                                    femmt_simulation_transformer_interleaved_core_fixed_loss_angle,
                                    femmt_simulation_transformer_integrated_core_fixed_loss_angle
                                    ]

        result_log_filepath_direction = os.path.join(os.path.dirname(__file__), "fixtures", "results")

        # femmt_simulation_inductor_core_material_database = os.path.join(result_log_filepath_direction, "log_electro_magnetic_inductor_core_material.json")
        # femmt_simulation_inductor_core_fixed_loss_angle = os.path.join(result_log_filepath_direction, "log_electro_magnetic_inductor_core_fixed_loss_angle.json")
        # femmt_simulation_inductor_core_fixed_loss_angle_litz_wire = os.path.join(result_log_filepath_direction, "log_electro_magnetic_inductor_core_fixed_loss_angle_litz_wire.json")
        # femmt_simulation_transformer_core_fixed_loss_angle = os.path.join(result_log_filepath_direction, "log_electro_magnetic_transformer_core_fixed_loss_angle.json")
        # femmt_simulation_transformer_interleaved_core_fixed_loss_angle = os.path.join(result_log_filepath_direction, "log_electro_magnetic_transformer_interleaved_core_fixed_loss_angle.json")
        # femmt_simulation_transformer_integrated_core_fixed_loss_angle = os.path.join(result_log_filepath_direction, "log_electro_magnetic_transformer_integrated_core_fixed_loss_angle.json")
        #
        #
        # fixture_log_filepath_list = []

        for count, filepath in enumerate(result_log_filepath_list):
            test_result_log = fmt.MagneticComponent.decode_settings_from_log(filepath, working_directory)

            assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

            #compare_result_logs(test_result_log, fixture_result_log)

def test_thermal_simulation(thermal_simulation):
    electromagnetoquasistatic_result_log, thermal_result_log = thermal_simulation

    # check magnetoquasistatic simulation results
    assert os.path.exists(electromagnetoquasistatic_result_log), "Electro magnetic simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "log_electro_magnetic_thermal.json")
    compare_result_logs(electromagnetoquasistatic_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "results_thermal.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)