import femmt as fmt
import numpy as np
import materialdatabase as mdb
import os

database = mdb.MaterialDatabase()
# ######################################   {RELUCTANCE_CALCULATION}   ##################################################
# Set core-geometry from core database
core_db = fmt.core_database()
# core_names = ["PQ 40/40", "PQ 20/20", "PQ 26/20", "PQ 26/25"]
core_names = ["PQ 40/40", "PQ 20/20"]

core_w_list = [core_db[name]["core_w"] for name in core_names]
window_h_list = [core_db[name]["window_h"] for name in core_names]
window_w_list = [core_db[name]["window_w"] for name in core_names]

# Set air-gap and core parameters
no_of_turns = list(np.linspace(8, 12, 5))                          # Set No. of turns (N)
n_air_gaps = list(np.linspace(1, 3, 3))                             # Set No. of air-gaps (n)
air_gap_height = list(np.linspace(0.0001, 0.0005, 2))            # Set air-gap length in metre (l)
air_gap_position = list(np.linspace(1, 99, 11))                     # Set air-gap position in percentage with respect to core window height

material_names = ["N95", "N87"]                                     # Set relative permeability in F/m (u)
mu_rel = [database.get_initial_permeability(material_name=name) for name in material_names]
# mu_rel = [2700, 3000, 3100, 3200]

# Set two types of equally distributed air-gaps (used only for air-gaps more than 1):
# Type 1: Equally distributed air-gaps including corner air-gaps (eg: air-gaps-position = [0, 50, 100] for 3 air-gaps)
# Type 2: Equally distributed air-gaps excluding corner air-gaps (eg: air-gaps-position = [25, 50, 75] for 3 air-gaps)
mult_air_gap_type = [2]  # 'Type1 = with corner air-gaps; 'Type2' = without air-gaps; 'Type0' = single air-gap

# Call to Reluctance model (Class MagneticCircuit)
mc1 = fmt.MagneticCircuit(core_w=core_w_list, window_h=window_h_list, window_w=window_w_list, no_of_turns=no_of_turns,
                          n_air_gaps=n_air_gaps, air_gap_h=air_gap_height, air_gap_position=air_gap_position,
                          mu_rel=mu_rel, mult_air_gap_type=mult_air_gap_type)

# Calculate total reluctance and creates data_matrix to access all corresponding parameters and results
mc1.core_reluctance()
mc1.air_gap_reluctance()
param = mc1.get_param_pos_dict()
n_cases_0 = len(mc1.data_matrix)
# To access all/any data from MagneticCircuit class, use self.data_matrix[:, parameter_position].
# The parameters are arranged as shown below:
# self.data_matrix = [core_w{0}, window_h{1}, window_w{2}, mu_rel{3}, no_of_turns{4}, n_air_gaps{5}, air_gap_h{6},
# air_gap_position{}7, mult_air_gap_type{8}, inductance{9}, core_h_middle{10}, r_inner{11}, r_outer{12},
# center_leg_area{13}, outer_leg_area{14}]
# Example: If you want to access inductance, type self.data_matrix[:, 9]

# ############################################   {FILTRATION}   #######################################################
# Applied constraints
goal_inductance = 120 * 1e-6
L_tolerance_percent = 10
conductor_radius = 0.0013
winding_factor = 0.91

current_max = 1                     # Max current amplitude with an assumption of sinusoidal current waveform
B_sat_N95 = 525 * 1e-3              # B_saturation for N95 material (in Tesla)
percent_of_B_sat = 90               # Percent of B_sat allowed in the designed core
percent_of_total_loss = 100          # Percent of total_loss allowed in FEM simulation
freq = 100 * 1e3                    # Switching frequency
mu_imag = 100
Cu_sigma = 5.96 * 1e7

# 1st Filter: Filter out cases where physical geometry is not possible
data_matrix_1 = mc1.data_matrix[np.where((mc1.data_matrix[:, param["no_of_turns"]] * np.pi * conductor_radius ** 2)
                < (winding_factor * mc1.data_matrix[:, param["window_h"]] * mc1.data_matrix[:, param["window_w"]]))]
n_cases_1 = len(data_matrix_1)

# 2nd Filter: Based on +-10% goal inductance tolerance band
data_matrix_2 = data_matrix_1[np.where((data_matrix_1[:, param["inductance"]] > ((100 - L_tolerance_percent) / 100) * goal_inductance) &
                                           (data_matrix_1[:, param["inductance"]] < ((100 + L_tolerance_percent) / 100) * goal_inductance))]
n_cases_2 = len(data_matrix_2)

total_flux_max = (data_matrix_2[:, param["inductance"]] * current_max) / data_matrix_2[:, param["no_of_turns"]]    # flux_max = L * i_max / N
data_matrix_2 = np.hstack((data_matrix_2, np.reshape(total_flux_max, (len(total_flux_max), 1))))  # position: 15
param["total_flux_max"] = 15

B_max_center = total_flux_max / data_matrix_2[:, param["center_leg_area"]]
B_max_middle = total_flux_max / (np.pi * data_matrix_2[:, param["core_w"]] * data_matrix_2[:, param["core_h_middle"]])
B_max_outer = total_flux_max / data_matrix_2[:, param["outer_leg_area"]]

data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_center, (len(B_max_center), 1))))  # position: 16
param["B_max_center"] = 16
data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_outer, (len(B_max_outer), 1))))  # position: 17
param["B_max_outer"] = 17

# 3rd Filter: Filter out cases where B_max is greater than B_sat
data_matrix_3 = data_matrix_2[np.where((B_max_center < (percent_of_B_sat / 100) * B_sat_N95) &
                                       (B_max_middle < (percent_of_B_sat / 100) * B_sat_N95) &
                                       (B_max_outer < (percent_of_B_sat / 100) * B_sat_N95))]
n_cases_3 = len(data_matrix_3)

# 4th Filter: Filter out data-matrix according to calculated hysteresis loss + DC winding loss
# Volume chosen as per "Masterthesis_Till_Piepenbrock" pg-45
volume_center = (np.pi * (data_matrix_3[:, param["core_w"]] / 2) ** 2) * \
                (data_matrix_3[:, param["window_h"]] + data_matrix_3[:, param["core_h_middle"]] -
                 (data_matrix_3[:, param["n_air_gaps"]] * data_matrix_3[:, param["air_gap_h"]]))
volume_outer = (np.pi * ((data_matrix_3[:, param["r_outer"]] ** 2) - (data_matrix_3[:, param["r_inner"]] ** 2))) * \
               (data_matrix_3[:, param["window_h"]] + data_matrix_3[:, param["core_h_middle"]])

P_hyst_center = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * ((data_matrix_3[:, param["B_max_center"]] /
                                                                 (fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2)
P_hyst_outer = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * ((data_matrix_3[:, param["B_max_outer"]] /
                                                                (fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2)

P_hyst_density_center = P_hyst_center * volume_center
P_hyst_density_middle = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * \
                        ((data_matrix_3[:, param["total_flux_max"]] / (fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2) * \
                        (1 / (2 * np.pi * data_matrix_3[:, param["core_h_middle"]])) * \
                        np.log((data_matrix_3[:, param["r_inner"]] * 2) / data_matrix_3[:, param["core_w"]])
P_hyst_density_outer = P_hyst_outer * volume_outer

total_hyst_loss = P_hyst_density_center + (2 * P_hyst_density_middle) + P_hyst_density_outer
data_matrix_3 = np.hstack((data_matrix_3, np.reshape(total_hyst_loss, (len(total_hyst_loss), 1)))) # position: 18
param["P_hyst_density_total"] = 18

# Winding loss (only DC loss)
Resistance = (data_matrix_3[:, param["no_of_turns"]] * 2 * np.pi *
              (data_matrix_3[:, param["core_w"]] / 2 + conductor_radius)) / \
             (Cu_sigma * (np.pi * (conductor_radius ** 2)))

DC_loss = ((current_max ** 2) / 2) * Resistance
data_matrix_3 = np.hstack((data_matrix_3, np.reshape(DC_loss, (len(DC_loss), 1))))  # position: 19
param["DC_loss"] = 19

total_loss = DC_loss + total_hyst_loss
data_matrix_3 = np.hstack((data_matrix_3, np.reshape(total_loss, (len(total_loss), 1))))  # position: 20
param["total_loss"] = 20

# Sort the data_matrix with respect to total losses column
data_matrix_3 = data_matrix_3[data_matrix_3[:, param["total_loss"]].argsort()]

FEM_data_matrix = data_matrix_3[0:int((percent_of_total_loss / 100) * len(data_matrix_3)), :]
n_cases_FEM = len(FEM_data_matrix)

# ##########################################   {FEM_SIMULATION}   ######################################################
qwerty = 1
FEM_inductance = np.zeros((len(FEM_data_matrix), 1))

# Component type
component = "inductor"

# Working directory can be set arbitrarily
working_directory = os.path.join(os.path.dirname(__file__), '..')

# MagneticComponent class object
geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory)

for i in range(len(FEM_data_matrix)):
    print(f"###################################{i}############################################")
    core = fmt.Core(core_w=FEM_data_matrix[i, param["core_w"]], window_w=FEM_data_matrix[i, param["window_w"]],
                    window_h=FEM_data_matrix[i, param["window_h"]],
                    material="95_100")
    # mu_rel=3000, phi_mu_deg=10,
    # sigma=0.5)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    if int(FEM_data_matrix[i, param["n_air_gaps"]]) == 1:
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_position"]],
                             FEM_data_matrix[i, param["air_gap_h"]])
    else:
        if int(FEM_data_matrix[i, param["mult_air_gap_type"]]) == 1:
            position_list = list(np.linspace(0, 100, int(FEM_data_matrix[i, param["n_air_gaps"]])))
            for position in position_list:
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, position,
                                     FEM_data_matrix[i, param["air_gap_h"]])

        elif int(FEM_data_matrix[i, param["mult_air_gap_type"]]) == 2:
            position_list = list(np.linspace(0, 100, int(FEM_data_matrix[i, param["n_air_gaps"]]) + 2))
            position_list.remove(0.0)
            position_list.remove(100.0)
            for position in position_list:
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, position,
                                     FEM_data_matrix[i, param["air_gap_h"]])

    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters: use solid wires
    winding = fmt.Winding(int(FEM_data_matrix[i, param["no_of_turns"]]), 0, fmt.Conductivity.Copper,
                          fmt.WindingType.Primary, fmt.WindingScheme.Square)
    winding.set_solid_conductor(0.0013)
    # winding.set_litz_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6, fill_factor=None)
    geo.set_windings([winding])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.001, 0.001)
    isolation.add_winding_isolations(0.0005)
    geo.set_isolation(isolation)

    # 5. create the model
    geo.create_model(freq=100000, visualize_before=False, save_png=False)

    # 6. start simulation
    geo.single_simulation(freq=100000, current=[3], show_results=False)

    FEM_inductance[i, 0] = geo.read_log()["single_sweeps"][0]["winding1"]["self_inductance"][0]

print(FEM_inductance)
