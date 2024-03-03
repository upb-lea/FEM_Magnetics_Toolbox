"""Contains integration tests to benchmark new changes against previous simulation results."""
import pytest
import os
import numpy as np
import json
import femmt as fmt
import deepdiff
import femmt.examples.basic_inductor
import femmt.examples.basic_transformer_interleaved
import femmt.examples.basic_transformer
import femmt.examples.basic_transformer_three_winding
import femmt.examples.basic_transformer_integrated
import femmt.examples.basic_transformer_center_tapped
import femmt.examples.basic_transformer_stacked
import femmt.examples.basic_transformer_stacked_center_tapped
import femmt.examples.basic_inductor_foil_vertical
import femmt.examples.basic_transformer_n_winding
import femmt.examples.advanced_inductor_sweep
import femmt.examples.basic_transformer_5_windings
import femmt.examples.basic_inductor_time_domain
import femmt.examples.basic_transformer_time_domain
import femmt.examples.advanced_inductor_air_gap_sweep
import femmt.examples.component_study.transformer_component_study
import materialdatabase as mdb


def compare_result_logs(first_log_filepath, second_log_filepath, significant_digits=6):
    """Compare the result log against a given one to see the differences when running the integration tests."""
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

    difference = deepdiff.DeepDiff(first_content, second_content, ignore_order=True,
                                   significant_digits=significant_digits)
    print(f"{difference = }")

    assert not deepdiff.DeepDiff(first_content, second_content, ignore_order=True,
                                 significant_digits=significant_digits)
    # made several tests with the deepdiff command:
    # tried adding not existing keys in one of the dicts: results as expected in an error
    # changed values in very nested dict: results as expected in an error
    # So this command is valid to compare the dicts.


def compare_thermal_result_logs(first_log_filepath, second_log_filepath, significant_digits=6):
    """Compare the thermal result log against a given one to see the differences when running the integration tests."""
    with open(first_log_filepath, "r") as fd:
        first_content = json.loads(fd.read())

    with open(second_log_filepath, "r") as fd:
        second_content = json.loads(fd.read())

    difference = deepdiff.DeepDiff(first_content, second_content, ignore_order=True,
                                   significant_digits=significant_digits)
    print(f"{difference = }")

    assert not deepdiff.DeepDiff(first_content, second_content, ignore_order=True,
                                 significant_digits=significant_digits)
    # made several tests with the deepdiff command:
    # tried adding not existing keys in one of the dicts: results as expected in an error
    # changed values in very nested dict: results as expected in an error
    # So this command is valid to compare the dicts.


@pytest.fixture
def temp_folder():
    """Fixture to create the temporary folder the results are stored."""
    # Setup temp folder
    temp_folder_path = os.path.join(os.path.dirname(__file__), "temp")

    if not os.path.exists(temp_folder_path):
        os.mkdir(temp_folder_path)

    config_path = os.path.join(os.path.dirname(__file__), "..", "config.json")
    if os.path.exists(config_path):
        with open(config_path, "r") as fd:
            content = json.load(fd)
            onelab_path = content["onelab"]
    elif os.path.isdir(os.path.join(os.path.dirname(__file__), "..", "..", "onelab")):
        onelab_path = os.path.join(os.path.dirname(__file__), "..", "..", "onelab")
    else:
        print("FEMMT Simulations will not work without onelab path, which is not set yet. In order to \
              run the tests please first specify a onelab path by creating a config.json in the tests folder.")
        onelab_path = None

    # Test
    yield temp_folder_path, onelab_path


@pytest.fixture
def femmt_simulation_inductor_core_material_database(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder
    
    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions, material=mdb.Material.N95, temperature=25, frequency=100000,
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
        winding.set_solid_round_conductor(conductor_radius=0.0013,
                                          conductor_arrangement=fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 9, None)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=100000, pre_visualize_geometry=False, save_png=False)

        geo.single_simulation(freq=100000, current=[4.5], show_fem_simulation_results=False)

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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme, colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_inductor_core_material_database_measurement(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions, material=mdb.Material.N95, temperature=25, frequency=100000,
                        permeability_datasource=fmt.MaterialDataSource.Measurement,
                        permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                        permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                        permittivity_datasource=fmt.MaterialDataSource.Measurement,
                        permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                        permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme, colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        mu_r_abs=3000, phi_mu_deg=10, sigma=0.5, permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme, colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle_litz_wire(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3000, phi_mu_deg=10, sigma=0.5, permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme,
                               colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical(temp_folder):
    """Fixture for the integration test."""
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
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

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

        vww.set_winding(winding, 5, fmt.WindingScheme.FoilVertical, wrap_para_type)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme,
                               colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal(temp_folder):
    """Fixture for the integration test."""
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
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme,
                               colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_transformer_core_fixed_loss_angle(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015,
                                                        window_w=0.012,
                                                        window_h=0.0295,
                                                        core_h=0.05)

        # 2. set core parameters
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=1.2, permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme,
                               colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_transformer_interleaved_core_fixed_loss_angle(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015,
                                                        window_w=0.012,
                                                        window_h=0.0295,
                                                        core_h=0.05)

        # 2. set core parameters
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        non_linear=False, sigma=1, re_mu_rel=3200, phi_mu_deg=10,
                        permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme,
                               colors_geometry=colors_geometry)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_transformer_integrated_core_fixed_loss_angle(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                    working_directory=working_directory, verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.02,
                                                        window_w=0.011,
                                                        window_h=0.03,
                                                        core_h=0.05)

        # 2. set core parameters
        core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                        mu_r_abs=3100, phi_mu_deg=12,
                        sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
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
                               case_gap_bot, show_thermal_simulation_results=False, pre_visualize_geometry=False,
                               color_scheme=color_scheme,
                               colors_geometry=colors_geometry,
                               flag_insulation=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_transformer_stacked_center_tapped(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer,
                                    working_directory=working_directory, verbosity=fmt.Verbosity.Silent,
                                    is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_dimensions = fmt.dtos.StackedCoreDimensions(core_inner_diameter=0.02, window_w=0.015, window_h_top=0.005,
                                                         window_h_bot=0.017)
        core = fmt.Core(core_type=fmt.CoreType.Stacked, core_dimensions=core_dimensions, mu_r_abs=3100, phi_mu_deg=12,
                        sigma=1.2,
                        permeability_datasource=fmt.MaterialDataSource.Custom,
                        permittivity_datasource=fmt.MaterialDataSource.Custom)
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Stacked, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.002, stacked_position=fmt.StackedPosition.Top)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, stacked_position=fmt.StackedPosition.Bot)
        geo.set_air_gaps(air_gaps)

        # set_center_tapped_windings() automatically places the condu
        insulation, coil_window, transformer_window = fmt.functions_topologies.set_center_tapped_windings(
            core=core,
            primary_turns=14,
            primary_radius=1.1e-3,
            primary_number_strands=50,
            primary_strand_radius=0.00011,
            secondary_parallel_turns=2,
            secondary_thickness_foil=1e-3,
            iso_top_core=0.001,
            iso_bot_core=0.001,
            iso_left_core=0.002,
            iso_right_core=0.001,
            iso_primary_to_primary=2e-4,
            iso_secondary_to_secondary=2e-4,
            iso_primary_to_secondary=4e-4,
            interleaving_type=fmt.CenterTappedInterleavingType.TypeC,
            interleaving_scheme=fmt.InterleavingSchemesFoilLitz.ter_3_4_sec_ter_4_3_sec,
            primary_coil_turns=3,
            primary_additional_bobbin=1e-3,
            winding_temperature=100,
            bobbin_coil_left=3e-3,
            center_foil_additional_bobbin=0e-3)

        geo.set_insulation(insulation)
        geo.set_winding_windows([coil_window, transformer_window])

        geo.create_model(freq=200000, pre_visualize_geometry=False)

        geo.single_simulation(freq=200000, current=[20, 120, 120], phi_deg=[0, 180, 180],
                              show_fem_simulation_results=False)

        geo.get_inductances(I0=1, op_frequency=200000)

        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given
        # magnetic component
        # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a
        # boundary temperature
        # which will be applied on the boundary of the simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": {  # epoxy resign
                "top": 1.54,
                "top_right": 1.54,
                "right": 1.54,
                "bot_right": 1.54,
                "bot": 1.54
            },
            "core": 5,  # ferrite
            "winding": 400,  # copper
            "air_gaps": 180,  # aluminium nitride
            "insulation": 0.42  # polyethylene
        }

        # Here the case size can be determined
        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
        # This does not change the results of the simulation (at least when every boundary is set equally) but will set
        # the temperature offset.
        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

        # Here the boundary sides can be turned on (1) or off (0)
        # By turning off the flag a neumann boundary will be applied at this point with heat flux = 0
        boundary_flags = {
            "flag_boundary_top": 0,
            "flag_boundary_top_right": 0,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        # In order for the thermal simulation to work an electro_magnetic simulation has to run before.
        # The em-simulation will create a file containing the losses.
        # When the losses file is already created and contains the losses for the current model, it is
        # enough to run geo.create_model in
        # order for the thermal simulation to work (geo.single_simulation is not needed).
        # Obviously when the model is modified and the losses can be out of date and therefore the
        # geo.single_simulation needs to run again.
        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, show_thermal_simulation_results=False,
                               color_scheme=fmt.colors_ba_jonas,
                               colors_geometry=fmt.colors_geometry_ba_jonas,
                               flag_insulation=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result


@pytest.fixture
def femmt_simulation_transformer_5_windings(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=16.1e-3, window_w=(22.5 - 12) / 2 * 1e-3,
                                                        core_inner_diameter=12e-3, core_h=22e-3)
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
        cells = winding_window.NHorizontalAndVerticalSplit(horizontal_split_factors=[0.48, 0.75],
                                                           vertical_split_factors=[None, [0.5, 0.85], None])

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
        geo.create_model(freq=100000, pre_visualize_geometry=False)
        geo.single_simulation(freq=100000, current=[1.625, 6.9, 4.9, 0.5, 1],
                              phi_deg=[0, 180, 180, 90, 89], show_fem_simulation_results=False)

        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the
        # given magnetic component. In order to use the thermal simulation, thermal conductivities for each material
        # can be entered as well as a boundary temperature which will be applied on the boundary
        # of the simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": {  # epoxy resign
                "top": 1.54,
                "top_right": 1.54,
                "right": 1.54,
                "bot_right": 1.54,
                "bot": 1.54
            },
            "core": 5,  # ferrite
            "winding": 400,  # copper
            "air_gaps": 180,  # aluminium nitride
            "insulation": 0.42  # polyethylene
        }

        # Here the case size can be determined
        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
        # This does not change the results of the simulation (at least when every boundary is set equally) but will
        # set the temperature offset.
        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

        # Here the boundary sides can be turned on (1) or off (0)
        # By turning off the flag a neumann boundary will be applied at this point with heat flux = 0
        boundary_flags = {
            "flag_boundary_top": 0,
            "flag_boundary_top_right": 0,
            "flag_boundary_right_top": 1,
            "flag_boundary_right": 1,
            "flag_boundary_right_bottom": 1,
            "flag_boundary_bottom_right": 1,
            "flag_boundary_bottom": 1
        }

        # In order for the thermal simulation to work an electro_magnetic simulation has to run before.
        # The em-simulation will create a file containing the losses.
        # When the losses file is already created and contains the losses for the current model, it is enough to run
        # geo.create_model in order for the thermal simulation to work (geo.single_simulation is not needed).
        # Obviously when the model is modified and the losses can be out of date and therefore the geo.single_simulation
        # needs to run again.
        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, show_thermal_simulation_results=False,
                               color_scheme=fmt.colors_ba_jonas,
                               colors_geometry=fmt.colors_geometry_ba_jonas,
                               flag_insulation=True)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    thermal_result = os.path.join(temp_folder_path, "results", "results_thermal.json")

    return electromagnetoquasistatic_result, thermal_result

@pytest.fixture
def femmt_simulation_inductor_time_domain(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Inductor,
                                    working_directory=working_directory, verbosity=fmt.Verbosity.Silent, is_gui=True)
        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])
        inductor_frequency = 270000
        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        material=mdb.Material.N49, temperature=45, frequency=inductor_frequency,
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
        winding.parallel = False

        # 7. add conductor to vww and add winding window to MagneticComponent
        vww.set_winding(winding, 7, None)
        geo.set_winding_windows([winding_window])

        # 8. create the model
        geo.create_model(freq=270000, pre_visualize_geometry=False, save_png=False)

        # 6.a. start simulation
        # time value list
        t = np.linspace(0, 1 / inductor_frequency, 5)
        t_list = [float(x) for x in t.tolist()]
        #  Current values list
        current_values = 4.5 * np.cos(2 * np.pi * inductor_frequency * t)
        current_values_list = current_values.tolist()

        geo.time_domain_simulation(current_period_vec=[current_values_list],
                                   time_period_vec=t_list,
                                   number_of_periods=1,
                                   show_fem_simulation_results=False,
                                   show_rolling_average=False,)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

    return electromagnetoquasistatic_result

@pytest.fixture
def femmt_simulation_transformer_time_domain(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Transformer,
                                    working_directory=working_directory, verbosity=fmt.Verbosity.Silent, is_gui=True)
        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015, window_w=0.012, window_h=0.0295, core_h=0.04)
        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        material=mdb.Material.N49, temperature=45, frequency=200000,
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

        winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
        winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)
        winding2.parallel = False

        # 7. add conductor to vww and add winding window to MagneticComponent
        bot.set_winding(winding2, 10, None)
        top.set_winding(winding1, 10, None)
        geo.set_winding_windows([winding_window])

        # 8. create the model
        geo.create_model(freq=200000, pre_visualize_geometry=False)

        # 6.a. start simulation
        # time value list
        t = np.linspace(0, 2 / 200000, 5)
        t_list = [float(x) for x in t.tolist()]
        # Current values
        current_values_1 = 2 * np.cos(2 * np.pi * 200000 * t)
        current_values_2 = 2 * np.cos(2 * np.pi * 200000 * t + np.pi)
        current_values_list_1 = current_values_1.tolist()
        current_values_list_2 = current_values_2.tolist()

        geo.time_domain_simulation(current_period_vec=[current_values_list_1, current_values_list_2],
                                   time_period_vec=t_list,
                                   number_of_periods=1,
                                   show_fem_simulation_results=False,
                                   show_rolling_average=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

    return electromagnetoquasistatic_result
@pytest.fixture
def femmt_simulation_transformer_3_windings_time_domain(temp_folder):
    """Fixture for the integration test."""
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Transformer,
                                    working_directory=working_directory, verbosity=fmt.Verbosity.Silent, is_gui=True)
        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_dimensions = fmt.dtos.SingleCoreDimensions(window_h=0.06, window_w=0.03, core_inner_diameter=0.015, core_h=0.08)
        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        material=mdb.Material.N49, temperature=45, frequency=200000,
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
        winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
        winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
        winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        winding3 = fmt.Conductor(2, fmt.Conductivity.Copper)
        winding3.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        # 7. add conductor to vww and add winding window to MagneticComponent
        top_left.set_winding(winding1, 7, fmt.WindingType.Single)
        top_right.set_winding(winding2, 8, fmt.WindingType.Single)
        bot_right.set_winding(winding3, 9, fmt.WindingType.Single)
        geo.set_winding_windows([winding_window])

        # 8. start simulation with given frequency, currents and phases
        geo.create_model(freq=250000, pre_visualize_geometry=False)

        # 6.a. start simulation
        # time value list
        t = np.linspace(0, 2 / 250000, 5)
        t_list = [float(x) for x in t.tolist()]
        # # Current values
        current_values_1 = 4 * np.cos(2 * np.pi * 250000 * t)
        current_values_2 = 4 * np.cos(2 * np.pi * 250000 * t + np.pi)
        current_values_3 = 4 * np.cos(2 * np.pi * 250000 * t)
        current_values_list_1 = current_values_1.tolist()
        current_values_list_2 = current_values_2.tolist()
        current_values_list_3 = current_values_3.tolist()

        geo.time_domain_simulation(current_period_vec=[current_values_list_1, current_values_list_2, current_values_list_3],
                                   time_period_vec=t_list,
                                   number_of_periods=1,
                                   show_fem_simulation_results=False,
                                   show_rolling_average=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")

    return electromagnetoquasistatic_result


def test_inductor_core_material_database(femmt_simulation_inductor_core_material_database):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_inductor_core_material_database

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_material.json")
    print(test_result_log)
    print(fixture_result_log)
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_inductor_core_material_database.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_inductor_core_material_database_measurement(femmt_simulation_inductor_core_material_database_measurement):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_inductor_core_material_database_measurement

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_material_measurement.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_inductor_core_material_database_measurement.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_inductor_core_fixed_loss_angle(femmt_simulation_inductor_core_fixed_loss_angle):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_inductor_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_inductor_core_fixed_loss_angle.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_inductor_core_fixed_loss_angle_litz_wire(femmt_simulation_inductor_core_fixed_loss_angle_litz_wire):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_inductor_core_fixed_loss_angle_litz_wire

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_fixed_loss_angle_litz_wire.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=3)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_inductor_core_fixed_loss_angle_litz_wire.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log, significant_digits=3)


def test_inductor_core_fixed_loss_angle_foil_vertical(femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_fixed_loss_angle_foil_vertical.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results

    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_inductor_core_fixed_loss_angle_foil_vertical.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_inductor_core_fixed_loss_angle_foil_horizontal(
        femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_inductor_core_fixed_loss_angle_foil_horizontal.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_inductor_core_fixed_loss_angle_foil_horizontal.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_transformer_core_fixed_loss_angle(femmt_simulation_transformer_core_fixed_loss_angle):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_transformer_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_transformer_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_transformer_core_fixed_loss_angle.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_transformer_interleaved_core_fixed_loss_angle(femmt_simulation_transformer_interleaved_core_fixed_loss_angle):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_transformer_interleaved_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_transformer_interleaved_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_transformer_interleaved_core_fixed_loss_angle.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_transformer_integrated_core_fixed_loss_angle(femmt_simulation_transformer_integrated_core_fixed_loss_angle):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_transformer_integrated_core_fixed_loss_angle

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "log_electro_magnetic_transformer_integrated_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_transformer_integrated_core_fixed_loss_angle.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)


def test_simulation_transformer_stacked_center_tapped(femmt_simulation_transformer_stacked_center_tapped):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_transformer_stacked_center_tapped

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "transformer_stacked_center_tapped.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "thermal_transformer_stacked_center_tapped.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log, significant_digits=3)


def test_simulation_transformer_5_windings(femmt_simulation_transformer_5_windings):
    """Integration test to validate the magnetoquasistatic simulation and the thermal simulation."""
    test_result_log, thermal_result_log = femmt_simulation_transformer_5_windings

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "transformer_5_windings.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)

    # check thermal simulation results
    # assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    # fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
    #                                   "thermal_transformer_5_windings.json")
    # compare_thermal_result_logs(thermal_result_log, fixture_result_log, significant_digits=2)

def test_simulation_inductor_time_domain(femmt_simulation_inductor_time_domain):
    """Integration test to validate the magnetoquasistatic simulation."""
    test_result_log = femmt_simulation_inductor_time_domain

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "inductor_time_domain.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)


def test_simulation_transformer_time_domain(femmt_simulation_transformer_time_domain):
    """Integration test to validate the magnetoquasistatic simulation."""
    test_result_log = femmt_simulation_transformer_time_domain

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "transformer_time_domain.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)


def test_simulation_transformer_3_windings_time_domain(femmt_simulation_transformer_3_windings_time_domain):
    """Integration test to validate the magnetoquasistatic simulation."""
    test_result_log = femmt_simulation_transformer_3_windings_time_domain

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results",
                                      "transformer_3_windings_time_domain.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)


def test_load_files(temp_folder, femmt_simulation_inductor_core_material_database,
                    femmt_simulation_inductor_core_fixed_loss_angle,
                    femmt_simulation_inductor_core_fixed_loss_angle_litz_wire,
                    femmt_simulation_inductor_core_fixed_loss_angle_foil_horizontal,
                    femmt_simulation_inductor_core_fixed_loss_angle_foil_vertical,
                    femmt_simulation_transformer_core_fixed_loss_angle,
                    femmt_simulation_transformer_interleaved_core_fixed_loss_angle,
                    femmt_simulation_transformer_integrated_core_fixed_loss_angle
                    ):
    """Function tests if simulations can be set up from a simulation file.

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

        # femmt_simulation_inductor_core_material_database = os.path.join(result_log_filepath_direction,
        # "log_electro_magnetic_inductor_core_material.json")
        # femmt_simulation_inductor_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "log_electro_magnetic_inductor_core_fixed_loss_angle.json")
        # femmt_simulation_inductor_core_fixed_loss_angle_litz_wire = os.path.join(result_log_filepath_direction,
        # "log_electro_magnetic_inductor_core_fixed_loss_angle_litz_wire.json")
        # femmt_simulation_transformer_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "log_electro_magnetic_transformer_core_fixed_loss_angle.json")
        # femmt_simulation_transformer_interleaved_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "log_electro_magnetic_transformer_interleaved_core_fixed_loss_angle.json")
        # femmt_simulation_transformer_integrated_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "log_electro_magnetic_transformer_integrated_core_fixed_loss_angle.json")
        #
        #
        # fixture_log_filepath_list = []

        for filepath in result_log_filepath_list:
            test_result_log = fmt.MagneticComponent.decode_settings_from_log(filepath, working_directory)

            assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

            # compare_result_logs(test_result_log, fixture_result_log)

##############################
# Basic example tests
# These tests just run the basic examples an see if the run without error
# There is no result comparison
##############################


def test_basic_examples(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor.basic_example_inductor(onelab_folder=onelab_folder,
                                                         show_visual_outputs=False,
                                                         is_test=True)


def test_basic_example_transformer_interleaved(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_interleaved.basic_example_transformer_interleaved(onelab_folder=onelab_folder,
                                                                                       show_visual_outputs=False,
                                                                                       is_test=True)


def test_basic_example_transformer(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer.basic_example_transformer(onelab_folder=onelab_folder, show_visual_outputs=False,
                                                               is_test=True)


def test_basic_example_transformer_three_winding(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_three_winding.basic_example_transformer_three_winding(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_example_transformer_integrated(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_integrated.basic_example_transformer_integrated(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_example_transformer_center_tapped(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_center_tapped.basic_example_transformer_center_tapped(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_example_transformer_stacked(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_stacked.basic_example_transformer_stacked(onelab_folder=onelab_folder,
                                                                               show_visual_outputs=False,
                                                                               is_test=True)


def test_basic_example_transformer_stacked_center_tapped(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_stacked_center_tapped.basic_example_transformer_stacked_center_tapped(
        onelab_folder=onelab_folder, show_visual_outputs=False, is_test=True)


def test_basic_example_inductor_foil_vertical(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor_foil_vertical.basic_example_inductor_foil_vertical(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_example_transformer_n_winding(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_n_winding.basic_example_transformer_n_winding(onelab_folder=onelab_folder,
                                                                                   show_visual_outputs=False,
                                                                                   is_test=True)


def test_basic_example_transformer_5_windings(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_5_windings.basic_example_transformer_5_windings(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)

def test_basic_inductor_time_domain(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor_time_domain.basic_example_inductor_time_domain(onelab_folder=onelab_folder,
                                                                                 show_visual_outputs=False,
                                                                                 is_test=True)

def test_basic_transformer_time_domain(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_time_domain.basic_example_transformer_time_domain(onelab_folder=onelab_folder,
                                                                                       show_visual_outputs=False,
                                                                                       is_test=True)

def test_basic_transformer_3_windings_time_domain(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_three_winding.basic_example_transformer_three_winding(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)
def test_advanced_example_inductor_sweep(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.advanced_inductor_sweep.advanced_example_inductor_sweep(onelab_folder=onelab_folder,
                                                                           show_visual_outputs=False,
                                                                           is_test=True)


def test_advanced_example_inductor_air_gap_sweep(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.advanced_inductor_air_gap_sweep.basic_example_sweep(onelab_folder=onelab_folder,
                                                                       show_visual_outputs=False,
                                                                       is_test=True)

def test_transformer_component_study(temp_folder):
    """Integrationtest to test the basic example file."""
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.component_study.transformer_component_study.transformer_component_study(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)
