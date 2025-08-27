"""Contains integration tests to benchmark new changes against previous simulation results.

This file is split up in two file:
 - part 1 contains only a few rudimentary tests: Run this to get a rudimentary, fast overview if your changes are ok
 - part 2 contains all the other tests: Run this to get a better overview if your changes are ok

To get the full picture, you need to run all tests!
"""

import json
import os

# 3rd party libraries
import deepdiff
import pytest
import numpy as np

# own libraries
import femmt as fmt


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


def compare_result_logs(first_log_filepath: str, second_log_filepath: str, significant_digits: int = 6,
                        ignore_order: bool = True, number_format_notation: str = "f"):
    """
    Compare the result log against a given one to see the differences when running the integration tests.

    :param first_log_filepath: filepath to the first log file
    :type first_log_filepath: str
    :param second_log_filepath: filepath to the second log file
    :type second_log_filepath: str
    :param significant_digits: significant digits to compare
    :type significant_digits: int
    :param ignore_order: True to ignore the order of the dict keys. Defaults to True.
    :type ignore_order: bool
    :param number_format_notation: number format notation for deepdiff
    :type number_format_notation: str
    """
    first_content = None
    second_content = None

    with open(first_log_filepath, "r") as fd:
        first_content = json.loads(fd.read())
        if "simulation_settings" in first_content:
            if "date" in first_content["simulation_settings"]:
                del first_content["simulation_settings"]["date"]
            if "working_directory" in first_content["simulation_settings"]:
                del first_content["simulation_settings"]["working_directory"]

    with open(second_log_filepath, "r") as fd:
        second_content = json.loads(fd.read())
        if "simulation_settings" in second_content:
            if "date" in second_content["simulation_settings"]:
                del second_content["simulation_settings"]["date"]
            if "working_directory" in second_content["simulation_settings"]:
                del second_content["simulation_settings"]["working_directory"]

    difference = deepdiff.DeepDiff(first_content, second_content, ignore_order=ignore_order,
                                   significant_digits=significant_digits, number_format_notation=number_format_notation)
    print(f"{difference=}")

    assert not deepdiff.DeepDiff(first_content, second_content, ignore_order=ignore_order,
                                 significant_digits=significant_digits, number_format_notation=number_format_notation)
    # made several tests with the deepdiff command:
    # tried adding not existing keys in one of the dicts: results as expected in an error
    # changed values in very nested dict: results as expected in an error
    # So this command is valid to compare the dicts.


def compare_thermal_result_logs(first_log_filepath: str, second_log_filepath: str, significant_digits: int = 6):
    """
    Compare the thermal result log against a given one to see the differences when running the integration tests.

    :param first_log_filepath: filepath to the first log file
    :type first_log_filepath: str
    :param second_log_filepath: filepath to the second log file
    :type second_log_filepath: str
    :param significant_digits: significant digits to compare
    :type significant_digits: int
    """
    with open(first_log_filepath, "r") as fd:
        first_content = json.loads(fd.read())

    with open(second_log_filepath, "r") as fd:
        second_content = json.loads(fd.read())

    difference = deepdiff.DeepDiff(first_content, second_content, ignore_order=True,
                                   significant_digits=significant_digits)
    print(f"{difference=}")

    assert not deepdiff.DeepDiff(first_content, second_content, ignore_order=True,
                                 significant_digits=significant_digits)
    # made several tests with the deepdiff command:
    # tried adding not existing keys in one of the dicts: results as expected in an error
    # changed values in very nested dict: results as expected in an error
    # So this command is valid to compare the dicts.

@pytest.fixture
def fixture_inductor_core_material_database_measurement(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    onelab_verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

        core_material = fmt.ImportedComplexCoreMaterial(material=fmt.Material.N49,
                                                        temperature=25,
                                                        permeability_datasource=fmt.DataSource.TDK_MDT,
                                                        permittivity_datasource=fmt.DataSource.LEA_MTB)

        core = fmt.Core(material=core_material,
                        core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        detailed_core_model=False)

        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([[0.0005]], per_layer_of_turns=False)
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=25)
        winding.set_solid_round_conductor(conductor_radius=0.0013,
                                          conductor_arrangement=fmt.ConductorArrangement.Square)

        vww.set_winding(winding, 9, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
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
    geometry_result = os.path.join(temp_folder_path, "results", "log_coordinates_description.json")
    material_result = os.path.join(temp_folder_path, "results", "log_material.json")

    return electromagnetoquasistatic_result, thermal_result, geometry_result, material_result

@pytest.fixture
def fixture_transformer_core_fixed_loss_angle(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    onelab_verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=0.015,
                                                        window_w=0.012,
                                                        window_h=0.0295,
                                                        core_h=0.05)

        core_material = fmt.LinearComplexCoreMaterial(mu_r_abs=3100,
                                                      phi_mu_deg=12,
                                                      dc_conductivity=1.2,
                                                      eps_r_abs=0,
                                                      phi_eps_deg=0)

        core = fmt.Core(material=core_material,
                        core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        detailed_core_model=False)

        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
        geo.set_air_gaps(air_gaps)

        # 4. set insulation
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
        insulation.add_winding_insulations([[0.0002, 0.0002], [0.0002, 0.0002]], per_layer_of_turns=False)
        geo.set_insulation(insulation)

        # 5. create winding window and virtual winding windows (vww)
        winding_window = fmt.WindingWindow(core, insulation)
        left, right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit, 0.0005)

        # 6. create conductors and set parameters
        winding1 = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=25)
        winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        winding2 = fmt.Conductor(1, fmt.ConductorMaterial.Copper, temperature=25)
        winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

        # 7. add conductor to vww and add winding window to MagneticComponent
        left.set_winding(winding1, 10, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
        right.set_winding(winding2, 10, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
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
    geometry_result = os.path.join(temp_folder_path, "results", "log_coordinates_description.json")
    material_result = os.path.join(temp_folder_path, "results", "log_material.json")

    return electromagnetoquasistatic_result, thermal_result, geometry_result, material_result

@pytest.fixture
def fixture_inductor_time_domain(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.TimeDomain, component_type=fmt.ComponentType.Inductor,
                                    working_directory=working_directory, onelab_verbosity=fmt.Verbosity.Silent, is_gui=True)
        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])

        core_material = fmt.LinearComplexCoreMaterial(mu_r_abs=3000,
                                                      phi_mu_deg=0,
                                                      dc_conductivity=1,
                                                      eps_r_abs=0,
                                                      phi_eps_deg=0)

        core = fmt.Core(material=core_material,
                        core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        detailed_core_model=False)

        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
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
        winding.parallel = False

        # 7. add conductor to vww and add winding window to MagneticComponent
        vww.set_winding(winding, 7, None, fmt.Align.ToEdges, fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
        geo.set_winding_windows([winding_window])

        # 8. create the model
        inductor_frequency = 270000
        geo.create_model(freq=inductor_frequency, pre_visualize_geometry=False, save_png=False)

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
                                   show_rolling_average=False, )

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electromagnetoquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_magnetic.json")
    geometry_result = os.path.join(temp_folder_path, "results", "log_coordinates_description.json")
    material_result = os.path.join(temp_folder_path, "results", "log_material.json")

    return electromagnetoquasistatic_result, geometry_result, material_result


@pytest.fixture
def fixture_inductor_electrostatic(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.ElectroStatic, component_type=fmt.ComponentType.Inductor,
                                    working_directory=working_directory, onelab_verbosity=fmt.Verbosity.Silent, is_gui=True)
        # Set onelab path manually
        geo.file_data.onelab_folder_path = onelab_folder

        # 2. set core parameters
        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                        window_w=core_db["window_w"],
                                                        window_h=core_db["window_h"],
                                                        core_h=core_db["core_h"])
        core_material = fmt.ElectrostaticCoreMaterial(eps_r=100e3)

        core = fmt.Core(material=core_material,
                        core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        detailed_core_model=False)

        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 50)
        # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0002, 90)
        geo.set_air_gaps(air_gaps)

        # 4. set insulation and the materials of insulation (optional)
        bobbin_db = fmt.bobbin_database()["PQ 40/40"]
        bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
                                                      bobbin_window_w=bobbin_db["bobbin_window_w"],
                                                      bobbin_window_h=bobbin_db["bobbin_window_h"],
                                                      bobbin_h=bobbin_db["bobbin_h"])
        insulation = fmt.Insulation(flag_insulation=True, bobbin_dimensions=bobbin_dimensions)
        bobbin_material = fmt.insulation_materials_database()["core_insulation"]["bobbins"]["Thermoset"]["Phenolic"]
        insulation.add_core_insulations(0.2e-3, 0.2e-3, 0.2e-3, 0.2e-3,
                                        dielectric_constant=bobbin_material["dielectric_constant"])
        turn_insulation_material = fmt.insulation_materials_database()["wire_insulation"]["plastic_insulation"]["Plenum Polyvinyl Chloride (Plenum PVC)"]
        insulation.add_turn_insulation([0.2e-3], dielectric_constant=[turn_insulation_material["dielectric_constant"]])
        # This is an air between turns if needed
        insulation.add_winding_insulations([[1e-3, 1e-3]], per_layer_of_turns=True)
        # Kapton material is added between every layer of turns
        layer_insulation = fmt.insulation_materials_database()["film_insulation"]["Kapton"]
        insulation.add_insulation_between_layers(thickness=0.6e-3, dielectric_constant=layer_insulation["dielectric_constant"])
        geo.set_insulation(insulation)

        # 5. create winding window and virtual winding windows (vww)
        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        # 6. create conductor and set parameters: use solid wires
        winding = fmt.Conductor(0, fmt.ConductorMaterial.Copper, temperature=45)
        winding.set_solid_round_conductor(conductor_radius=1.1506e-3, conductor_arrangement=fmt.ConductorArrangement.Square)
        winding.parallel = False  # set True to make the windings parallel, currently only for solid conductors
        # 7. add conductor to vww and add winding window to MagneticComponent
        vww.set_winding(winding, 7, None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.VerticalUpward_HorizontalRightward,
                        zigzag=True)
        geo.set_winding_windows([winding_window])

        # 8. create the model
        geo.create_model(freq=2700000, pre_visualize_geometry=False, save_png=False, skin_mesh_factor=0.5)

        # 8. run electrostatic simulation
        num_turns_w1 = 7
        V_A = 1
        V_B = 0
        voltages_winding_1 = [
            V_A - (V_A - V_B) * i / (num_turns_w1 - 1)
            for i in range(num_turns_w1)
        ]
        geo.electrostatic_simulation(voltage=[voltages_winding_1], ground_outer_boundary=False, core_voltage=0,
                                     show_fem_simulation_results=False, save_to_excel_file=False)

    except Exception as e:
        print("An error occurred while creating the femmt mesh files:", e)
    except KeyboardInterrupt:
        print("Keyboard interrupt..")

    electrostaticquasistatic_result = os.path.join(temp_folder_path, "results", "log_electro_static.json")
    geometry_result = os.path.join(temp_folder_path, "results", "log_coordinates_description.json")
    material_result = os.path.join(temp_folder_path, "results", "log_material.json")

    return electrostaticquasistatic_result, geometry_result, material_result

@pytest.fixture
def fixture_planar_transformer_interleaved(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder

    # Create new temp folder, build model and simulate
    try:
        working_directory = temp_folder_path
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory,
                                    onelab_verbosity=fmt.Verbosity.Silent, is_gui=True)

        # Set onelab path manually
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
        insulation.add_winding_insulations([[1.15e-4, 1.15e-4], [1.15e-4, 1.15e-4]])
        insulation.add_insulation_between_layers(thickness=0.01e-3)
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        winding1 = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
        winding1.set_rectangular_conductor(thickness=35e-6, width=1.524e-3)

        winding2 = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
        winding2.set_rectangular_conductor(thickness=35e-6, width=4.8e-3)
        winding2.parallel = True

        vww.set_interleaved_winding(winding1, 5, winding2, 2, fmt.InterleavedWindingScheme.VerticalAlternating,
                                    foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalUpward, group_size=1)
        geo.set_winding_windows([winding_window])

        geo.create_model(freq=1000000, pre_visualize_geometry=False)

        geo.single_simulation(freq=1000000, current=[3, 5], phi_deg=[0, 180], show_fem_simulation_results=False)

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
    geometry_result = os.path.join(temp_folder_path, "results", "log_coordinates_description.json")

    return electromagnetoquasistatic_result, thermal_result, geometry_result


def test_inductor_core_material_database_measurement(fixture_inductor_core_material_database_measurement: pytest.fixture):
    """
    Integration test to validate the magnetoquasistatic simulation and the thermal simulation.

    :param fixture_inductor_core_material_database_measurement: fixture for inductor
    :type fixture_inductor_core_material_database_measurement: pytest.fixture
    """
    test_result_log, thermal_result_log, geometry_result_log, material_result_log = fixture_inductor_core_material_database_measurement

    assert os.path.exists(material_result_log), "Material log creation did not work!"

    fixture_material_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                        "material_inductor_core_material_database_measurement.json")
    compare_result_logs(material_result_log, fixture_material_log, significant_digits=2, ignore_order=False, number_format_notation="e")

    assert os.path.exists(geometry_result_log), "Geometry creation did not work!"

    fixture_geometry_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                        "geometry_inductor_core_material_measurement.json")
    compare_result_logs(geometry_result_log, fixture_geometry_log, significant_digits=6)

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "log_inductor_core_material_measurement.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=5)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "thermal_inductor_core_material_database_measurement.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log, significant_digits=6)

def test_transformer_core_fixed_loss_angle(fixture_transformer_core_fixed_loss_angle: pytest.fixture):
    """
    Integration test to validate the magnetoquasistatic simulation and the thermal simulation.

    :param fixture_transformer_core_fixed_loss_angle: fixture for the transformer
    :type fixture_transformer_core_fixed_loss_angle: pytest.fixture
    """
    test_result_log, thermal_result_log, geometry_result_log, material_result_log = fixture_transformer_core_fixed_loss_angle

    assert os.path.exists(geometry_result_log), "Geometry creation did not work!"

    fixture_geometry_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                        "geometry_transformer_core_fixed_loss_angle.json")
    compare_result_logs(geometry_result_log, fixture_geometry_log, significant_digits=10)

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "log_transformer_core_fixed_loss_angle.json")
    compare_result_logs(test_result_log, fixture_result_log)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "thermal_transformer_core_fixed_loss_angle.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)

def test_simulation_inductor_time_domain(fixture_inductor_time_domain: pytest.fixture):
    """
    Integration test to validate the magnetoquasistatic simulation.

    :param fixture_inductor_time_domain: fixture for the inductor time domain simulation
    :type fixture_inductor_time_domain: pytest.fixture
    """
    test_result_log, geometry_result_log, material_result_log = fixture_inductor_time_domain

    assert os.path.exists(geometry_result_log), "Geometry creation did not work!"

    fixture_geometry_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                        "geometry_inductor_time_domain.json")
    compare_result_logs(geometry_result_log, fixture_geometry_log, significant_digits=10)

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "log_inductor_time_domain.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)

def test_simulation_planar_transformer_interleaved(fixture_planar_transformer_interleaved: pytest.fixture):
    """
    Integration test to validate the magnetoquasistatic simulation.

    :param fixture_planar_transformer_interleaved: fixture for the inductor time domain simulation
    :type fixture_planar_transformer_interleaved: pytest.fixture
    """
    test_result_log, thermal_result_log, geometry_result_log = fixture_planar_transformer_interleaved

    assert os.path.exists(geometry_result_log), "Geometry creation did not work!"

    fixture_geometry_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                        "geometry_planar_transformer_interleaved.json")
    compare_result_logs(geometry_result_log, fixture_geometry_log, significant_digits=10)

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "log_planar_transformer_interleaved.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)

    # check thermal simulation results
    assert os.path.exists(thermal_result_log), "Thermal simulation did not work!"
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "thermal_planar_transformer_interleaved.json")
    compare_thermal_result_logs(thermal_result_log, fixture_result_log)

def test_simulation_inductor_electrostatic(fixture_inductor_electrostatic: pytest.fixture):
    """
    Integration test to validate the electrostatic simulation.

    :param fixture_inductor_electrostatic: fixture for the inductor electrostatic  simulation
    :type fixture_inductor_electrostatic: pytest.fixture
    """
    test_result_log, geometry_result_log, material_result_log = fixture_inductor_electrostatic

    assert os.path.exists(geometry_result_log), "Geometry creation did not work!"

    fixture_geometry_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                        "geometry_inductor_electrostatic.json")
    compare_result_logs(geometry_result_log, fixture_geometry_log, significant_digits=10)

    assert os.path.exists(test_result_log), "Electro static simulation did not work!"

    # e_s mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures",
                                      "inductor_electrostatic.json")
    compare_result_logs(test_result_log, fixture_result_log, significant_digits=4)

def test_load_files(temp_folder: pytest.fixture,
                    fixture_inductor_core_material_database_measurement: pytest.fixture,
                    fixture_transformer_core_fixed_loss_angle: pytest.fixture,
                    fixture_inductor_time_domain: pytest.fixture,
                    fixture_inductor_electrostatic: pytest.fixture
                    ):
    """Function tests if simulations can be set up from a simulation file.

    There is no complete function check, there is just an error-check if the load will fail or not.

    Note: Fixtures are used, to make sure that the latest self-generated log-files can be read.
    It is not okay to use predefined logfiles which where generated once.

    :param temp_folder: fixture for temporary simulation folder and onelab folder
    :type temp_folder: pytest.fixture
    :param fixture_inductor_core_material_database_measurement: fixture for inductor
    :type fixture_inductor_core_material_database_measurement: pytest.fixture,
    :param fixture_transformer_core_fixed_loss_angle: fixture for transformer
    :type fixture_transformer_core_fixed_loss_angle: pytest.fixture
    :param fixture_inductor_time_domain: fixture
    :type fixture_inductor_time_domain: pytest.fixture
    :param fixture_inductor_electrostatic: fixture
    :type fixture_inductor_electrostatic: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder

    working_directory = temp_folder_path
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

        result_log_filepath_list = [fixture_inductor_core_material_database_measurement,
                                    fixture_transformer_core_fixed_loss_angle,
                                    fixture_inductor_time_domain,
                                    fixture_inductor_electrostatic
                                    ]

        result_log_filepath_direction = os.path.join(os.path.dirname(__file__), "fixtures", "results")

        # femmt_simulation_inductor_core_material_database = os.path.join(result_log_filepath_direction,
        # "inductor_core_material.json")
        # femmt_simulation_inductor_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "inductor_core_fixed_loss_angle.json")
        # femmt_simulation_inductor_core_fixed_loss_angle_litz_wire = os.path.join(result_log_filepath_direction,
        # "inductor_core_fixed_loss_angle_litz_wire.json")
        # femmt_simulation_transformer_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "transformer_core_fixed_loss_angle.json")
        # femmt_simulation_transformer_interleaved_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "transformer_interleaved_core_fixed_loss_angle.json")
        # femmt_simulation_transformer_integrated_core_fixed_loss_angle = os.path.join(result_log_filepath_direction,
        # "transformer_integrated_core_fixed_loss_angle.json")
        #
        #
        # fixture_log_filepath_list = []

        for filepath in result_log_filepath_list:
            test_result_log = fmt.MagneticComponent.decode_settings_from_log(filepath, working_directory)

            assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

            # compare_result_logs(test_result_log, fixture_result_log)
