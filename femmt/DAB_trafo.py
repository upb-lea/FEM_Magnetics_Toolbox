from FEMMT import MagneticComponent
import numpy as np


# -- Definitions --
inductance_calculation = 0
skin_accuracy = 0.5


# -- Parameters --
frequency = 200000
N1 = 22  # Turns in main window
N2 = 8  # Turns in main window
Ns1 = 12  # Turns in stray window
Ns2 = 0  # Turns in stray window


# -- FFT --
I0 = 4


# -- Component Design --
geo = MagneticComponent(component_type="transformer")
geo.update_core(core_type="EI", window_h=0.028, window_w=0.007, core_w=0.02)

geo.update_air_gaps(method="percent",
                    n_air_gaps=5,
                    position_tag=[0, 0, 0, 0, 0],
                    air_gap_h=[0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
                    air_gap_position=[20, 35, 55, 75, 95],
                    stray_path=[[1, 2], geo.core_w/2+geo.window_w-0.0005])

geo.update_conductors(n_turns=[[N1, Ns1], [N2, Ns2]],
                      conductor_type=["litz", "litz"],
                      winding=["interleaved", "interleaved"],
                      scheme=["horizontal", "horizontal"],
                      conductor_radix=[0.0007, 0.001],
                      litz_para_type=['implicite_strands_number', 'implicite_strands_number'],
                      ff=[0.5, 0.5],
                      strand_radix=[0.01e-3, 0.01e-3],
                      cond_cond_isolation=[0.0001, 0.0001, 0.0003],
                      core_cond_isolation=[0.0005])

# -- Simulation --
if inductance_calculation:
    geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)
else:
    geo.single_simulation(freq=frequency, current=[I0, N1/N2*I0], phi=[0, 180], skin_mesh_factor=skin_accuracy)
