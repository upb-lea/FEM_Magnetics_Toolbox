import pytest
import femmt
import numpy as np

def test_fft():
    example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28], [-175.69, 103.47, 175.69, -103.47, -175.69]])

    out = femmt.fft(example_waveform, plot='no', mode='rad', f0=25000, title='FFT input current')
    out_f_test = np.array([25013., 75039., 125065., 175091., 275143., ])
    out_x_test = np.array([169.6, 26.8, 3., 5.4, 2.])
    out_phi_test = np.array([-2.5, 2.8, -2.5, -3.1, -2.7])
    out_test = np.array([out_f_test, out_x_test, out_phi_test])
    assert out.any() == out_test.any()

    out = femmt.fft(example_waveform, plot='no', mode='rad', f0=25000, title='FFT input current', filter_type='harmonic', filter_value_harmonic=10)
    out_f_test = np.array([0., 25013., 50026.])
    out_x_test = np.array([0.2, 169.6, 0.2])
    out_phi_test = np.array([3.1, -2.5, 1.8])
    out_test = np.array([out_f_test, out_x_test, out_phi_test])
    assert out.any() == out_test.any()

    with pytest.raises(ValueError):
        femmt.fft(example_waveform, mode='deg')
        femmt.fft(example_waveform, mode='rad')
        femmt.fft(example_waveform, mode='unallowed_mode')

def test_find_common_frequencies():
    frequency_1 = [50, 100, 150, 200]
    frequency_2 = [50, 100, 150, 170, 200]
    amplitude_1 = [1, 2, 3, 4]
    amplitude_2 = [5, 6, 7, 8, 9]
    phase_1 = [10, 20, 30, 40]
    phase_2 = [50, 60, 70, 80, 90]

    out = femmt.find_common_frequencies(frequency_1, amplitude_1, phase_1, frequency_2, amplitude_2, phase_2)
    common_f = [50, 100, 150, 200]
    common_a = [[1, 5], [2, 6], [3, 7], [4, 9]]
    common_phase = [[10, 50], [20, 60], [30, 70], [40, 90]]
    out_test = [common_f, common_a, common_phase]
    assert out == out_test

def test_cost_function_core():
    assert femmt.cost_function_core(1, "ferrite") == 5.5
    assert femmt.cost_function_core(1.25, "ferrite") == 6.875
    assert femmt.cost_function_core(1.25, "nanocristalline") == 28.75


def test_phases_deg_from_time_current():
    time_vec = [0, 1, 2, 3, 4]
    current_1 = [0, 1, 3, 2, 1]
    current_2 = [0, 2, 1, 3, 2]

    phase_deg_1, phase_deg_2 = femmt.phases_deg_from_time_current(time_vec, current_1, current_2)

    assert phase_deg_1 == 180
    assert phase_deg_2 == 270

def test_max_value_from_value_vec():
    current_1 = [0, 1, 3, 2, 1]
    current_2 = [0, 2, 1, 3, 2]

    max_1, max_2 = femmt.max_value_from_value_vec(current_1, current_2)

    assert max_1 == 3
    assert max_2 == 3


def test_cost_function_winding():
    assert femmt.cost_function_winding([0.1], [femmt.ConductorType.RoundSolid.name]) == [4.7]
    assert femmt.cost_function_winding([0.1, 0.9], [femmt.ConductorType.RoundSolid.name, femmt.ConductorType.RectangularSolid.name]) == [4.7, 35.1]
    single_strand_cross_section_75um = (75e-6 / 2) ** 2 * np.pi
    assert femmt.cost_function_winding([0.1], [femmt.ConductorType.RoundLitz.name], [single_strand_cross_section_75um]) == [pytest.approx(7, rel=1e-3)]


def test_cost_function_total():
    assert femmt.cost_function_total(1.25, "ferrite", [0.1, 0.9], [femmt.ConductorType.RoundSolid.name, femmt.ConductorType.RectangularSolid.name]) == \
           pytest.approx(62.233, rel=1e-3)


def test_reluctance():
    core_inner_diameter = 0.045
    single_air_gap_total_hight = 0.0002
    core_hight = 0.01
    tablet_hight = 0.01
    window_w = 0.02

    r_gap_round_round = femmt.r_air_gap_round_round(single_air_gap_total_hight, core_inner_diameter, core_hight, core_hight)
    assert r_gap_round_round == pytest.approx(97100, rel=1e-3)

    r_gap_round_inf = femmt.r_air_gap_round_inf(single_air_gap_total_hight, core_inner_diameter, core_hight)
    assert r_gap_round_inf == pytest.approx(94983, rel=1e-3)

    r_gap_tablet_cylinder = femmt.r_air_gap_tablet_cyl(tablet_hight, single_air_gap_total_hight, core_inner_diameter, window_w)
    assert r_gap_tablet_cylinder == pytest.approx(51694, rel=1e-3)

    r_gap_tablet_cylinder_real = femmt.r_air_gap_tablet_cyl_no_2d_axi(tablet_hight, single_air_gap_total_hight, core_inner_diameter, window_w)
    assert r_gap_tablet_cylinder_real == pytest.approx(82517, rel=1e-3)


def test_volume():
    window_h = 0.03
    window_w = 0.011
    core_inner_diameter = 0.02
    assert femmt.calculate_core_2daxi_total_volume(core_inner_diameter, window_h, window_w) == pytest.approx(6.798406502368311e-05, rel=1e-3)


def test_rms_current_calculation():
    time_current = [[0, 0.5e-6, 2.5e-6, 3e-6, 5e-6], [16.55, -10.55, -16.55, 10.55, 16.55]]
    assert femmt.i_rms(time_current) == pytest.approx(12.779756127041967, rel=1e-3)

def test_winding_resistance_calculation_solid():
    core_inner_diameter = 20e-3
    window_w = 0.05
    turns_count = 20
    conductor_radius = 1e-3
    material = "Copper"

    assert femmt.resistance_solid_wire(core_inner_diameter, window_w, turns_count, conductor_radius, material) == pytest.approx(0.02251034482758622, rel=1e-3)

def test_calculate_inductance_matrix():
    l_s_target_value = 85e-6
    l_h_target_value = 850e-6
    n_target_value = 3.1

    inductance_matrix = femmt.calculate_inductance_matrix_from_ls_lh_n(l_s_target_value, l_h_target_value, n_target_value)

    assert inductance_matrix[0] == pytest.approx([0.000935, 0.00027419354838709674], rel=1e-3)
    assert inductance_matrix[1] == pytest.approx([0.00027419354838709674, 8.844953173777314e-05], rel=1e-3)

def test_conductivity_temperature():
    temperature = 100

    copper_sigma_100_degree_calculated = femmt.conductivity_temperature("Copper", temperature)
    aluminium_sigma_100_degree_calculated = femmt.conductivity_temperature("Aluminium", temperature)

    assert copper_sigma_100_degree_calculated == pytest.approx(4.4874e7, rel=1e-3)
    assert aluminium_sigma_100_degree_calculated == pytest.approx(2.8627e7, rel=1e-3)
