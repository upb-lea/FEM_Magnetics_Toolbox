# minimal example for github clone
from femmt import MagneticComponent
import numpy as np

# component = "inductor"
component = "transformer"
# component = "integrated_transformer"


# Create Object
if component == "inductor":
    geo = MagneticComponent(component_type="inductor")

    # Update Geometry
    geo.core.update(type="EI", window_h=0.03, window_w=0.011)

    # geo.air_gaps.update(method="percent", n_air_gaps=4, air_gap_h=[0.0005, 0.0005, 0.0005, 0.0005],
    #                     position_tag=[0, 0, 0, 0], air_gap_position=[20, 40, 60, 80])
    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.0005], position_tag=[0])

    geo.update_conductors(n_turns=[[14]], conductor_type=["solid"], conductor_radii=[0.0015],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.0005], cond_cond_isolation=[0.0001])

    # geo.single_simulation(freq=100000, current=[1])
    geo.femm_reference(freq=100000, current=[1], sigma_cu=58, sign=[1], non_visualize=0)


if component == "transformer":
    simulate_before_thermal = True

    geo = MagneticComponent(component_type="transformer")
    geo.visualize_before = False

    # Update Geometry
    geo.core.update(type="EI", window_h=0.0295, window_w=0.012, core_w=0.015)

    # geo.air_gaps.update(n_air_gaps=0)
    geo.air_gaps.update(method="percent", n_air_gaps=1, air_gap_h=[0.0005],
                        air_gap_position=[50], position_tag=[0])

    geo.update_conductors(n_turns=[[36], [11]], conductor_type=["solid", "litz"],
                        litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                        ff=[None, 0.6], strands_numbers=[None, 600], strand_radii=[70e-6, 35.5e-6],
                        conductor_radii=[0.0011, None],
                        winding=["interleaved"], scheme=["horizontal"],
                        core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

    if simulate_before_thermal:
        # Perform a single simulation
        geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180,
                                                                                    170])
        # geo.single_simulation(freq=250000, current=[4.18368713, 4.28975166], phi_deg=[-1.09710805/np.pi*180,
        #                                                                               - 1.47917789/np.pi*180 + 180])

        # geo.get_inductances(I0=8, op_frequency=250000, skin_mesh_factor=0.5)
        # geo.femm_reference(freq=100000, current=[1, 2], sigma_cu=58, sign=[1, -1], non_visualize=0)

    thermal_conductivity_dict = {
            "air": 0.0263,
            "case": 0.3,
            "core": 5,
            "winding": 400 
    }

    geo.thermal_simulation(thermal_conductivity_dict)
    #geo.femm_thermal_validation(thermal_conductivity_dict)

if component == "integrated_transformer":
    geo = MagneticComponent(component_type="integrated_transformer")

    geo.update_conductors(n_turns=[[1, 3], [2, 6]], conductor_type=["litz", "litz"],
                          litz_para_type=['implicit_litz_radius', 'implicit_litz_radius'],
                          ff=[0.5, 0.5], strands_numbers=[100, 100], strand_radii=[70e-6, 70e-6],
                          winding=["interleaved", "interleaved"], scheme=["horizontal", "horizontal"],
                          core_cond_isolation=[0.0005, 0.0005], cond_cond_isolation=[0.0002, 0.0002, 0.0005])

    geo.core.update(type="EI", window_h=0.03, window_w=0.011)

    geo.air_gaps.update(method="percent",
                        n_air_gaps=2,
                        position_tag=[0, 0],
                        air_gap_h=[0.001, 0.001],
                        air_gap_position=[30, 40])

    geo.stray_path.update(start_index=0,
                          radius=geo.core.core_w / 2 + geo.core.window_w - 0.001)

    # Perform a single simulation
    geo.single_simulation(freq=250000, current=[8.0, 4.0], phi_deg=[0, 180])
    # geo.get_inductances(I0=10, op_frequency=100000, skin_mesh_factor=0.5)


# fft:
#
# [{'power': 2000, 'voltage': 1400, 'frequency': 250000, 'wp_lambda': 2.8,
# 'wp_n': 0.52, 'ratio_lm_ls': 10, 'wp_ib_il_vec': [[0,
# 0.07102806898251754, 3.141592653589793, 3.2126207225723107,
# 6.283185307179586], [-0.4232418074018094, 5.3261170692099515,
# 0.4232418074018094, -5.3261170692099515, -0.4232418074018094]],
# 'wp_ib_il_frequency_vec': array([  250000.,   750000.,  1250000.,
# 1750000.,  2250000.,  2750000.,
#           3250000.,  3750000.,  4250000.,  4750000.,  5250000., 5750000.,
#           6250000.,  6750000.,  7250000.,  7750000.,  8250000., 8750000.,
#           9250000.,  9750000., 10250000., 10750000.]),
# 'wp_ib_il_amplitude_vec': array([4.18368713, 1.23724602, 0.73192879,
# 0.51860308, 0.4001389 ,
#          0.32441466, 0.27161389, 0.23254221, 0.20234873, 0.17823164,
#          0.15845962, 0.14190533, 0.12780332, 0.11561575, 0.1049538 ,
#          0.09552949, 0.08712508, 0.07957295, 0.07274202, 0.06652836,
#          0.06084852, 0.05563475]), 'wp_ib_il_phase_rad_vec':
# array([-1.09710805, -1.48992872, -1.6308925 , -1.73079185, -1.81675942,
#          -1.8963576 , -1.97251841, -2.04661581, -2.11937918, -2.19123181,
#          -2.26243635, -2.33316464, -2.40353403, -2.47362759, -2.54350598,
#          -2.61321476, -2.68278905, -2.7522566 , -2.8216399 , -2.89095765,
#          -2.9602259 , -3.02945879])}

#   {'power': 2000, 'voltage': 1400, 'frequency': 250000, 'wp_lambda':
# 2.8, 'wp_n': 0.52, 'ratio_lm_ls': 10, 'wp_ib_il_vec': [[0,
# 0.07102806898251754, 3.141592653589793, 3.2126207225723107,
# 6.283185307179586], [-6.206083824802166, -6.5, 6.206083824802166, 6.5,
# -6.206083824802166]], 'wp_ib_ilm_frequency_vec': array([ 250000.,
# 750000., 1250000., 1750000.]), 'wp_ib_ilm_amplitude_vec':
# array([5.2733732 , 0.58539888, 0.21036846, 0.10705326]),
# 'wp_ib_ilm_phase_rad_vec': array([3.07377326, 2.93804   , 2.80202566,
# 2.66555066])}

#   {'power': 2000, 'voltage': 1400, 'frequency': 250000, 'wp_lambda':
# 2.8, 'wp_n': 0.52, 'ratio_lm_ls': 10, 'wp_ob_il_vec': [[0,
# 0.07102806898251754, 3.141592653589793, 3.2126207225723107,
# 6.283185307179586], [3.0070778490481853, 6.149580875989175,
# -3.0070778490481853, -6.149580875989175, 3.0070778490481853]],
# 'wp_ob_il_frequency_vec': array([ 250000.,  750000., 1250000., 1750000.,
# 2250000., 2750000.,
#          3250000., 3750000., 4250000., 4750000., 5250000., 5750000.,
#          6250000., 6750000.]), 'wp_ob_il_amplitude_vec':
# array([4.28975166, 0.78517217, 0.42402217, 0.29181263, 0.22233953,
#          0.17909771, 0.14938799, 0.12760058, 0.11086301, 0.09754825,
#          0.08666436, 0.0775716 , 0.06983873, 0.06316436]),
# 'wp_ob_il_phase_rad_vec': array([-0.51739939, -1.1086375 , -1.38030193,
# -1.54846645, -1.67539063,
#          -1.7822816 , -1.87804012, -1.96700668, -2.05153927, -2.13302653,
#          -2.212339  , -2.29005047, -2.36655512, -2.44213341])}
