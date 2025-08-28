"""Basic example testing.

These tests just run the basic examples and see if the run without error.
There is no result comparison.
"""

# python libraries
import os
import json

# 3rd party libraries
import pytest

# own libraries
import femmt.examples.basic_inductor
import femmt.examples.basic_transformer_interleaved
import femmt.examples.basic_transformer
import femmt.examples.basic_transformer_three_winding
import femmt.examples.basic_transformer_integrated
import femmt.examples.basic_transformer_center_tapped
import femmt.examples.basic_transformer_stacked
import femmt.examples.basic_transformer_stacked_center_tapped
import femmt.examples.basic_inductor_foil_vertical
import femmt.examples.basic_inductor_foil
import femmt.examples.basic_transformer_n_winding
import femmt.examples.advanced_inductor_sweep
import femmt.examples.basic_transformer_5_windings
import femmt.examples.basic_transformer_6_windings
import femmt.examples.experimental_inductor_time_domain
import femmt.examples.experimental_transformer_time_domain
import femmt.examples.experimental_transformer_three_winding_time_domain
import femmt.examples.basic_inductor_electrostatic
import femmt.examples.basic_transformer_electrostatic
import femmt.examples.advanced_inductor_air_gap_sweep
import femmt.examples.component_study.transformer_component_study
import femmt.examples.basic_transformer_excitation_sweep
import femmt.examples.basic_inductor_excitation_sweep
import femmt.examples.basic_split_windings
import femmt.examples.basic_planar_transformer
import femmt.examples.basic_planar_transformer_interleaved
import femmt.examples.PCB_planar.planar_transformer


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

def test_basic_example_inductor(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor.basic_example_inductor(onelab_folder=onelab_folder,
                                                         show_visual_outputs=False,
                                                         is_test=True)


def test_basic_example_transformer_interleaved(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_interleaved.basic_example_transformer_interleaved(onelab_folder=onelab_folder,
                                                                                              show_visual_outputs=False,
                                                                                              is_test=True)


def test_basic_example_transformer(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer.basic_example_transformer(onelab_folder=onelab_folder, show_visual_outputs=False,
                                                               is_test=True)


def test_basic_example_transformer_three_winding(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_three_winding.basic_example_transformer_three_winding(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_example_transformer_integrated(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_integrated.basic_example_transformer_integrated(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_example_transformer_center_tapped(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_center_tapped.basic_example_transformer_center_tapped(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_example_transformer_stacked(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_stacked.basic_example_transformer_stacked(onelab_folder=onelab_folder,
                                                                               show_visual_outputs=False,
                                                                               is_test=True)


def test_basic_example_transformer_stacked_center_tapped(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_stacked_center_tapped.basic_example_transformer_stacked_center_tapped(
        onelab_folder=onelab_folder, show_visual_outputs=False, is_test=True)


def test_basic_example_inductor_foil_vertical(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor_foil_vertical.basic_example_inductor_foil_vertical(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_example_inductor_foil(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor_foil.basic_example_inductor_foil(onelab_folder=onelab_folder,
                                                                   show_visual_outputs=False,
                                                                   is_test=True)


def test_basic_example_transformer_n_winding(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_n_winding.basic_example_transformer_n_winding(onelab_folder=onelab_folder,
                                                                                   show_visual_outputs=False,
                                                                                   is_test=True)


def test_basic_example_transformer_5_windings(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_5_windings.basic_example_transformer_5_windings(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_example_transformer_6_windings(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_6_windings.basic_example_transformer_6_windings(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_inductor_time_domain(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.experimental_inductor_time_domain.basic_example_inductor_time_domain(onelab_folder=onelab_folder,
                                                                                        show_visual_outputs=False,
                                                                                        is_test=True)


def test_basic_transformer_time_domain(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.experimental_transformer_time_domain.basic_example_transformer_time_domain(onelab_folder=onelab_folder,
                                                                                              show_visual_outputs=False,
                                                                                              is_test=True)


def test_basic_inductor_electrostatic(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor_electrostatic.basic_example_inductor_electrostatic(onelab_folder=onelab_folder,
                                                                                     show_visual_outputs=False,
                                                                                     is_test=True)


def test_basic_transformer_electrostatic(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_electrostatic.basic_example_transformer_electrostatic(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_transformer_3_windings_time_domain(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.experimental_transformer_three_winding_time_domain.basic_example_transformer_three_windings_time_domain(onelab_folder=onelab_folder,
                                                                                                                           show_visual_outputs=False,
                                                                                                                           is_test=True)


def test_advanced_example_inductor_sweep(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.advanced_inductor_sweep.advanced_example_inductor_sweep(onelab_folder=onelab_folder,
                                                                           show_visual_outputs=False,
                                                                           is_test=True)


def test_advanced_example_inductor_air_gap_sweep(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.advanced_inductor_air_gap_sweep.basic_example_sweep(onelab_folder=onelab_folder,
                                                                       show_visual_outputs=False,
                                                                       is_test=True)


def test_basic_transformer_component_study(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.component_study.transformer_component_study.transformer_component_study(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_inductor_excitation_sweep(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_inductor_excitation_sweep.basic_example_inductor_excitation_sweep(onelab_folder=onelab_folder,
                                                                                           show_visual_outputs=False,
                                                                                           is_test=True)


def test_basic_transformer_excitation_sweep(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_transformer_excitation_sweep.basic_example_transformer_excitation_sweep(onelab_folder=onelab_folder,
                                                                                                 show_visual_outputs=False,
                                                                                                 is_test=True)


def test_basic_transformer_split_windings(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    for num_windings in [2, 3, 5, 6]:
        print(f"Running simulation for {num_windings} windings")
        femmt.examples.basic_split_windings.run_transformer_vvw_split_examples(num_windings, onelab_folder=onelab_folder,
                                                                               show_visual_outputs=False, is_test=True)

def test_basic_planar_transformer(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_planar_transformer.basic_example_planar_transformer(onelab_folder=onelab_folder,
                                                                             show_visual_outputs=False,
                                                                             is_test=True)

def test_basic_planar_transformer_interleaved(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.basic_planar_transformer_interleaved.basic_example_planar_transformer_interleaved(onelab_folder=onelab_folder,
                                                                                                     show_visual_outputs=False,
                                                                                                     is_test=True)
def test_planar_transformer(temp_folder: pytest.fixture):
    """
    Integration test to the basic example file.

    :param temp_folder: temporary folder path and onelab filepath
    :type temp_folder: pytest.fixture
    """
    temp_folder_path, onelab_folder = temp_folder
    femmt.examples.PCB_planar.planar_transformer.basic_example_planar_transformer_using_vww(onelab_folder=onelab_folder,
                                                                                            show_visual_outputs=False,
                                                                                            is_test=True)
