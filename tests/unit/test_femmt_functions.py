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

    out = femmt.fft(example_waveform, plot='no', mode='rad', f0=25000, title='FFT input current', filter_type='harmonic', filter_value_harmonic = 10)
    out_f_test = np.array([0., 25013., 50026.])
    out_x_test = np.array([0.2, 169.6,   0.2])
    out_phi_test = np.array([3.1, -2.5,  1.8])
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

    out = femmt.find_common_frequencies(frequency_1, amplitude_1, phase_1, frequency_2,
                                                                   amplitude_2, phase_2)
    common_f = [50, 100, 150, 200]
    common_a = [[1, 5], [2, 6], [3, 7], [4, 9]]
    common_phase = [[10, 50], [20, 60], [30, 70], [40, 90]]
    out_test = [common_f, common_a, common_phase]
    assert out == out_test
