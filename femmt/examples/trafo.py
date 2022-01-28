from femmt import MagneticComponent
import numpy as np

# -- Definitions --
frequency = 200000
flag_reference = 0
inductance_calculation = 0
skin_accuracy = 0.5
N1 = 26  # Turns in main window
N2 = 24 # Turns in main window
Ns1 = 2  # Turns in stray window
Ns2 = 2  # Turns in stray window
I0 = 10

# -- Component Design --
geo = MagneticComponent(component_type="transformer")
geo.update_core(core_type="EI", window_h=0.02, window_w=0.01, core_w=0.02)
"""
geo.update_air_gaps(method="percent",
                    n_air_gaps=5,
                    position_tag=[0, 0, 0, 0, 0],
                    air_gap_h=[0.0005, 0.0005, 0.0005, 0.0005, 0.0005],
                    air_gap_position=[15, 35, 55, 75, 95])
"""
geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.002])
# geo.update_air_gaps(method="manually", n_air_gaps=0, air_gap_h=[0.002])
"""
geo.update_conductors(n_turns=[[N1, 0], [0, N2]],  # n_turns=[[N1, Ns1], [N2, Ns2]],
                      conductor_type=["litz", "solid"],
                      winding=["primary", "secondary"],
                      scheme=["hexa", "square"],
                      litz_para_type=['implicite_strands_number', 'implicite_strands_number'],
                      conductor_radix=[0.0010, 0.0008],
                      ff=[0.5, 0.5],
                      strand_radix=[0.01e-3, 0.01e-3],
                      thickness=[0.001, None],
                      wrap_para=["fixed_thickness", None],
                      cond_cond_isolation=[0.0001, 0.0001, 0.0003],
                      core_cond_isolation=[0.0005])
"""
geo.update_conductors(n_turns=[[N1], [N2]],  # n_turns=[[N1, Ns1], [N2, Ns2]],
                      conductor_type=["litz", "litz"],
                      winding=["interleaved"],
                      scheme=[["hexa", "hexa"]],
                      litz_para_type=['implicite_strands_number', 'implicite_strands_number'],
                      conductor_radix=[0.0010, 0.0008],
                      ff=[0.5, 0.5],
                      strand_radix=[0.01e-3, 0.01e-3],
                      thickness=[0.001, None],
                      wrap_para=["fixed_thickness", None],
                      cond_cond_isolation=[0.0001, 0.0001, 0.0003],
                      core_cond_isolation=[0.0005])

# -- Simulation --
if inductance_calculation and not flag_reference:
    geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)
else:
    if flag_reference:
        geo.femm_reference(freq=frequency, current=[I0, N1/N2*I0], sign=[1, -1], sigma=58, non_visualize=0)
        # geo.femm_reference(freq=frequency, current=[I0, 0], sign=[1, -1], sigma=58, non_visualize=0)
    else:
        geo.single_simulation(freq=frequency, current=[I0, N1/N2*I0], phi=[0, 180], skin_mesh_factor=skin_accuracy)

        # TODO: Warum reagiert get_Core_Loss() nicht mehr auf update_conductors() ?
        # geo.get_Core_Loss(Ipeak=[I0, N1/N2*I0], ki=0.53, alpha=1.50, beta=2.38, t_rise=0.3/frequency,  # according to
        #                  t_fall=0.3/frequency, f_switch=frequency, skin_mesh_factor=skin_accuracy)  # LEA Project 2020
