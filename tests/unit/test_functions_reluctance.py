"""Unit tests for the functions_reluctance.py module."""
# 3rd party libraries
import pytest

# own modules
import femmt.functions_reluctance as fr

def test_resistance_litz_wire():
    """Test function for resistance_litz_wire()."""
    # horizontal winding scheme
    resistance_horizontal = fr.resistance_litz_wire(
        0.01, window_w=0.01, window_h=0.012, turns_count=20, iso_core_top=1e-3, iso_core_bot=1e-3, iso_core_left=1e-3, iso_core_right=1e-3,
        iso_primary_to_primary=10e-6, temperature=100, scheme="horizontal_first", material="Copper", litz_wire_name="1.1x60x0.1")
    assert resistance_horizontal == pytest.approx(0.05772275287356322)

    # vertical winding scheme
    resistance_vertical = fr.resistance_litz_wire(
        0.01, window_w=0.01, window_h=0.012, turns_count=20, iso_core_top=1e-3, iso_core_bot=1e-3, iso_core_left=1e-3, iso_core_right=1e-3,
        iso_primary_to_primary=10e-6, temperature=100, scheme="vertical_first", material="Copper", litz_wire_name="1.1x60x0.1")
    assert resistance_vertical == pytest.approx(0.043211097701149434)
