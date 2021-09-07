
import FEMMT
import numpy as np
import unittest

class TestStringMethods(unittest.TestCase):

    def test_fft(self):
        example_waveform = np.array([[0, 1.34, 3.14, 4.48, 6.28], [-175.69, 103.47, 175.69, -103.47, -175.69]])
        out = FEMMT.fft(example_waveform, plot='no', rad='yes', f0=25000, title='FFT input current')
        out_f_test = np.array([25013., 75039., 125065., 175091., 275143., ])
        out_x_test = np.array([169.6, 26.8, 3., 5.4, 2.])
        out_phi_test = np.array([-2.5, 2.8, -2.5, -3.1, -2.7])
        out_test = np.array([out_f_test, out_x_test, out_phi_test])
        self.assertEqual(out.any(), out_test.any())

if __name__ == '__main__':
    unittest.main()