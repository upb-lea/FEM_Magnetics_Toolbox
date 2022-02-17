import pytest
import numpy as np
import os
import json
import shutil
from femmt import MagneticComponent

def compare_results(first_mesh, second_mesh):
    first_content = None
    second_content = None

    with open(first_mesh, "rb") as fd:
        first_content = fd.read(-1)

    with open(second_mesh, "rb") as fd:
        second_content = fd.read(-1)

    return first_content == second_content

@pytest.fixture
def temp_folder():
    # Setup
    temp_folder_path = os.path.join(os.path.dirname(__file__), "temp")

    if not os.path.exists(temp_folder_path):
        os.mkdir(temp_folder_path)

    # Create config.json
    config_file_path = os.path.join(os.path.dirname(__file__), "..", "..", "onelab")
    if os.path.exists(config_file_path):
        config_content = {"onelab": config_file_path}
        with open(os.path.join(temp_folder_path, "config.json"), "w") as fd:
            fd.write(json.dumps(config_content))
    else:
        assert 0, "onelab folder not found"
        exit(0)

    # Test
    yield temp_folder_path

    # Teardown
    # Removing the temp folder does not work yet.
    #shutil.rmtree(temp_folder_path)

@pytest.fixture
def femmt_simulation(temp_folder):
    # Create new temp folder, build model and simulate
    try:
        geo = MagneticComponent(component_type="transformer", working_directory=temp_folder)

        geo.core.update(window_h=0.0295, window_w=0.012, core_w=0.015)
        geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                            air_gap_position=[50], position_tag=[0])
        geo.update_conductors(n_turns=[[36], [11]], conductor_type=["solid", "litz"],
                            litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                            ff=[None, 0.6], strands_numbers=[None, 600], strand_radii=[70e-6, 35.5e-6],
                            conductor_radii=[0.0011, None],
                            winding=["interleaved"], scheme=["horizontal"],
                            core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

        geo.create_model(freq=250000, visualize_before=False)
        geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180, 170], show_results=False)

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

    return os.path.join(temp_folder, "results", "result_log_electro_magnetic.json")

def test_femmt(femmt_simulation):
    """
    The first idea was to compare the simulated meshes with test meshes simulated manually.
    It turns out that the meshes cannot be compared because even slightly differences in the mesh,
    can cause to a test failure, because the meshes are binary files.
    Those differences could even occur when running the simulation on different machines
    -> This was observed when creating a docker image and running the tests.

    Now as an exaple only the result log will be checked.
    """
    test_result_log = femmt_simulation

    assert os.path.exists(test_result_log), "Electro magnetic simulation did not work!"

    # e_m mesh
    fixture_result_log = os.path.join(os.path.dirname(__file__), "fixtures", "results", "result_log_electro_magnetic.json")
    assert compare_results(test_result_log, fixture_result_log), "Electro magnetic mesh is wrong."
