import femmt as fmt
import numpy as np
import re
import os
from os import listdir
from os.path import isfile, join
import shutil
import matplotlib.pyplot as plt
from itertools import product
import logging
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import proj3d
import materialdatabase as mdb
material_db = mdb.MaterialDatabase()


def automated_design_func():
    # ########################################   {DESIGN PARAMETERS}   #################################################
    save_directory_name = "sweep_new_3d_with_cost"  # New directory is created in FEM_Magnetics_Toolbox/femmt/examples/
    goal_inductance = 120 * 1e-6
    L_tolerance_percent = 10                # Inductance tolerance of +-10% is applied
    winding_factor = 0.91
    i_max = 8                               # Max current amplitude with assumption of sinusoidal current waveform 8
    percent_of_B_sat = 70                   # Percent of B_sat allowed in the designed core
    percent_of_total_loss = 30              # Percent of total_loss allowed in FEM simulation 30

    freq = 100 * 1e3                        # Switching frequency in Hz
    # mu_imag = 100                           # TODO: (completed)
    Cu_sigma = 5.96 * 1e7                   # copper conductivity (sigma) @ 20 degree celsius


    # Set core-geometry from core database or/and manual entry
    db_core_names = []  # "PQ 40/40", "PQ 40/30"    blank

    manual_core_w = list(np.linspace(0.005, 0.05, 10))# 10
    manual_window_h = list(np.linspace(0.01, 0.08, 5))# 5
    manual_window_w = list(np.linspace(0.005, 0.04, 10))# 10

    all_manual_combinations = list(product(manual_core_w, manual_window_h, manual_window_w))
    manual_core_w = [item[0] for item in all_manual_combinations]
    manual_window_h = [item[1] for item in all_manual_combinations]
    manual_window_w = [item[2] for item in all_manual_combinations]

    core_db = fmt.core_database()
    db_core_w = [core_db[core_name]["core_inner_diameter"] for core_name in db_core_names]
    db_window_h = [core_db[core_name]["window_h"] for core_name in db_core_names]
    db_window_w = [core_db[core_name]["window_w"] for core_name in db_core_names]

    core_w_list = db_core_w + manual_core_w
    window_h_list = db_window_h + manual_window_h
    window_w_list = db_window_w + manual_window_w

    # Set winding settings (Solid and Litz winding type)
    solid_conductor_r = [0.0013]  # blank
    litz_names = ["1.5x105x0.1"]  # "1.5x105x0.1", "1.4x200x0.071"

    litz_db = fmt.litz_database()
    litz_conductor_r = [litz_db[litz_name]["conductor_radii"] for litz_name in litz_names]
    litz_strand_r = [litz_db[litz_name]["strand_radii"] for litz_name in litz_names]
    litz_strand_n = [litz_db[litz_name]["strands_numbers"] for litz_name in litz_names]

    min_conductor_r = min(litz_conductor_r + solid_conductor_r)

    # Set air-gap and core parameters
    no_of_turns = [2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14,15,16,17,18,19,20]  # Set No. of turns (N)   2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14,15,16,17,18,19,20
    print(f"{no_of_turns = }")
    n_air_gaps = [1, 2]  # Set No. of air-gaps (n)
    air_gap_height = list(np.linspace(0.0001, 0.0005, 5))  # Set air-gap length in metre (l)
    air_gap_position = list(np.linspace(20, 80, 2))  # Set air-gap position in percent w.r.t. core window height

    material_names = ["N95"]  # Set relative permeability in F/m (u) , "N87"

    mu_rel = [material_db.get_material_property(material_name=material_name, property="initial_permeability")
              for material_name in material_names]

    temp_var = [material_db.permeability_data_to_pro_file(25, freq, material_name, "manufacturer_datasheet")
                for material_name in material_names]
    print(temp_var[0][1])
    mu_imag_interpol_func = [interp1d(temp_var[i][0], temp_var[i][1], kind="cubic") for i in range(len(material_names))]
    # print(f[0](0.2))
    # Set two types of equally distributed air-gaps (used only for air-gaps more than 1):
    # Type 1: Equally distributed air-gaps including corner air-gaps (eg: air-gaps-position = [0, 50, 100])
    # Type 2: Equally distributed air-gaps excluding corner air-gaps (eg: air-gaps-position = [25, 50, 75])
    # 'Type1 = with corner air-gaps; 'Type2' = without corner air-gaps; 'Type0' = single air-gap
    mult_air_gap_type = [2] #Type1-Edge, Type2: Centre
    # TODO: check if the issue has been resolved

    # ######################################   {RELUCTANCE_CALCULATION}   ##############################################
    # Call to Reluctance model (Class MagneticCircuit)
    mc1 = fmt.MagneticCircuit(core_w=core_w_list, window_h=window_h_list, window_w=window_w_list,
                              no_of_turns=no_of_turns, n_air_gaps=n_air_gaps, air_gap_h=air_gap_height,
                              air_gap_position=air_gap_position, mu_rel=mu_rel, mult_air_gap_type=mult_air_gap_type,
                              air_gap_method='percent', component_type='inductor', sim_type='sweep')
    param = mc1.get_param_pos_dict()
    n_cases_0 = len(mc1.data_matrix)
    print(n_cases_0)
    # Calculate total reluctance and creates data_matrix to access all corresponding parameters and results
    # To access all/any data from MagneticCircuit class, use self.data_matrix[:, param["parameter_name"]].
    # The parameters are arranged as shown below:
    # Example: If you want to access inductance, type self.data_matrix[:, param["inductance"]]

    # ############################################   {FILTRATION}   ####################################################
    # 1st Filter: ------------------------------------------------------------------------------------------------------
    # Filter out cases where physical geometry is not possible
    data_matrix_1 = mc1.data_matrix[np.where((mc1.data_matrix[:, param["no_of_turns"]] * np.pi * min_conductor_r ** 2)
                                             < (winding_factor * mc1.data_matrix[:,
                                                                 param["window_h"]] * mc1.data_matrix[:,
                                                                                      param[
                                                                                          "window_w"]]))]
    n_cases_1 = len(data_matrix_1)
    print(n_cases_1)
    # 2nd Filter:-------------------------------------------------------------------------------------------------------
    # Based on +-10% goal inductance tolerance band
    data_matrix_2 = data_matrix_1[
        np.where((data_matrix_1[:, param["inductance"]] > ((100 - L_tolerance_percent) / 100) * goal_inductance) &
                 (data_matrix_1[:, param["inductance"]] < ((100 + L_tolerance_percent) / 100) * goal_inductance))]
    n_cases_2 = len(data_matrix_2)
    print(n_cases_2)
    # 3rd Filter:-------------------------------------------------------------------------------------------------------
    # Filter out cases where B_max is greater than B_sat
    B_sat_dict = {}
    counter1 = 0
    for material_name in material_names:
        B_sat_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
        B_sat_dict[B_sat_key] = material_db.get_material_property(material_name=material_name, property="max_flux_density")
        counter1 = counter1 + 1

    # Creating B_saturated array corresponding to the material type
    B_sat = np.zeros((len(data_matrix_2), 1))
    for index in range(len(data_matrix_2)):
        B_sat[index] = B_sat_dict[data_matrix_2[index, param["mu_rel"]]]

    # flux_max = L * i_max / N
    total_flux_max = (data_matrix_2[:, param["inductance"]] * i_max) / data_matrix_2[:, param["no_of_turns"]]
    data_matrix_2 = np.hstack((data_matrix_2, np.reshape(total_flux_max, (len(total_flux_max), 1))))
    param["total_flux_max"] = 15

    B_max_center = total_flux_max / data_matrix_2[:, param["center_leg_area"]]
    B_max_middle = total_flux_max / (
            np.pi * data_matrix_2[:, param["core_w"]] * data_matrix_2[:, param["core_h_middle"]])
    B_max_outer = total_flux_max / data_matrix_2[:, param["outer_leg_area"]]

    data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_center, (len(B_max_center), 1))))
    param["B_max_center"] = 16
    data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_outer, (len(B_max_outer), 1))))
    param["B_max_outer"] = 17

    data_matrix_3 = np.zeros((0, 18))
    for index in range(len(data_matrix_2)):
        if (data_matrix_2[index, param["B_max_center"]] < (percent_of_B_sat / 100) * B_sat[index]) & \
                (data_matrix_2[index, param["B_max_outer"]] < (percent_of_B_sat / 100) * B_sat[index]):
            data_matrix_3 = np.vstack([data_matrix_3, data_matrix_2[index, :]])
    n_cases_3 = len(data_matrix_3)
    print(n_cases_3)
    # 4th Filter:-------------------------------------------------------------------------------------------------------
    # Filter out data-matrix according to calculated hysteresis loss + DC winding loss
    # Volume chosen as per "Masterthesis_Till_Piepenbrock" pg-45
    mu_imag_dict = {}
    counter_x = 0
    for material_name in material_names:
        mu_imag_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
        mu_imag_dict[mu_imag_key] = counter_x
        counter_x = counter_x + 1

    # Creating B_saturated array corresponding to the material type
    mu_imag = np.zeros(len(data_matrix_3))
    for index in range(len(data_matrix_3)):
        mu_imag[index] = mu_imag_interpol_func[mu_imag_dict[data_matrix_3[index, param["mu_rel"]]]](data_matrix_3[index, param["B_max_center"]])


    volume_center = (np.pi * (data_matrix_3[:, param["core_w"]] / 2) ** 2) * \
                    (data_matrix_3[:, param["window_h"]] + data_matrix_3[:, param["core_h_middle"]] -
                     (data_matrix_3[:, param["n_air_gaps"]] * data_matrix_3[:, param["air_gap_h"]]))
    volume_outer = (np.pi * ((data_matrix_3[:, param["r_outer"]] ** 2) - (data_matrix_3[:, param["r_inner"]] ** 2))) * \
                   (data_matrix_3[:, param["window_h"]] + data_matrix_3[:, param["core_h_middle"]])

    P_hyst_center = 0.5 * (2 * np.pi * freq) * fmt.mu0 * mu_imag * ((data_matrix_3[:, param["B_max_center"]] /
                                                                     (fmt.mu0 * data_matrix_3[:,
                                                                                param["mu_rel"]])) ** 2)
    # P_hyst_center = np.array(mu_imag) * np.array(P_hyst_center)
    P_hyst_outer = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * ((data_matrix_3[:, param["B_max_outer"]] /
                                                                    (fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2)

    P_hyst_density_center = P_hyst_center * volume_center
    P_hyst_density_middle = 0.5 * (2 * np.pi * freq) * mu_imag * fmt.mu0 * \
                            ((data_matrix_3[:, param["total_flux_max"]] / (
                                    fmt.mu0 * data_matrix_3[:, param["mu_rel"]])) ** 2) * \
                            (1 / (2 * np.pi * data_matrix_3[:, param["core_h_middle"]])) * \
                            np.log((data_matrix_3[:, param["r_inner"]] * 2) / data_matrix_3[:, param["core_w"]])
    P_hyst_density_outer = P_hyst_outer * volume_outer

    total_hyst_loss = P_hyst_density_center + (2 * P_hyst_density_middle) + P_hyst_density_outer
    data_matrix_3 = np.hstack((data_matrix_3, np.reshape(total_hyst_loss, (len(total_hyst_loss), 1))))  # position: 18
    param["P_hyst_density_total"] = 18

    # Winding loss (only DC loss)
    Resistance = (data_matrix_3[:, param["no_of_turns"]] * 2 * np.pi *
                  (data_matrix_3[:, param["core_w"]] / 2 + min_conductor_r)) / \
                 (Cu_sigma * (np.pi * (min_conductor_r ** 2)))

    DC_loss = ((i_max ** 2) / 2) * Resistance
    data_matrix_3 = np.hstack((data_matrix_3, np.reshape(DC_loss, (len(DC_loss), 1))))  # position: 19
    param["DC_loss"] = 19

    total_loss = DC_loss + total_hyst_loss
    data_matrix_3 = np.hstack((data_matrix_3, np.reshape(total_loss, (len(total_loss), 1))))  # position: 20
    param["total_loss"] = 20
    data_matrix_3 = data_matrix_3[data_matrix_3[:, param["total_loss"]].argsort()]
    data_matrix_4 = data_matrix_3[0:int((percent_of_total_loss / 100) * len(data_matrix_3)), :]

    total_loss = data_matrix_4[:, param["total_loss"]]
    max_total_loss = max(total_loss)
    normalized_total_loss = total_loss / max_total_loss
    data_matrix_4 = np.hstack((data_matrix_4, np.reshape(normalized_total_loss, (len(normalized_total_loss), 1))))  # position: 20
    param["normalized_total_loss"] = 21

    total_volume = np.pi * (data_matrix_4[:, param["r_outer"]] ** 2) * data_matrix_4[:, param["core_h_middle"]]
    max_volume = max(total_volume)
    normalized_total_volume = total_volume / max_volume
    data_matrix_4 = np.hstack((data_matrix_4, np.reshape(total_volume, (len(total_volume), 1))))  # position: 20
    param["total_volume"] = 22
    data_matrix_4 = np.hstack((data_matrix_4, np.reshape(normalized_total_volume, (len(normalized_total_volume), 1))))  # position: 20
    param["normalized_total_volume"] = 23

    # Sort the data_matrix with respect to total losses column----------------------------------------------------------
    # FEM_data_matrix = data_matrix_4[np.where(data_matrix_4[:, param["normalized_total_loss"]] <= ((0.01 / data_matrix_4[:, param["normalized_total_volume"]]) + 0.6))]
    FEM_data_matrix = data_matrix_4
    n_cases_FEM = len(FEM_data_matrix)
    print(n_cases_FEM)

    fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
    fmt.plt.title("Normalised volume vs Normalised losses")
    fmt.plt.xlabel("Normalised Volume")
    fmt.plt.ylabel("Normalised Losses")
    ax.plot(FEM_data_matrix[:, param['normalized_total_volume']], FEM_data_matrix[:, param['normalized_total_loss']], 'o')
    ax.grid()
    fmt.plt.show()

    # ##########################################   {FEM_SIMULATION}   ##################################################
    qwerty = 1

    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    working_directory = os.path.join(example_results_folder, "inductor")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    if not os.path.exists(os.path.join(os.path.dirname(__file__), save_directory_name)):
        os.mkdir(os.path.join(os.path.dirname(__file__), save_directory_name))

    working_directories = []
    file_names = []

    src_path = os.path.join(os.path.dirname(__file__), "example_results/inductor/results/log_electro_magnetic.json")

    counter3 = 0
    for j in range(len(solid_conductor_r) + len(litz_conductor_r)):
        conductor_r_list = litz_conductor_r + solid_conductor_r
        for i in range(len(FEM_data_matrix)):
            print(f"value of i:{i}")
            print(f"value of j:{j}")
            if not ((FEM_data_matrix[i, param["no_of_turns"]] * np.pi * conductor_r_list[j] ** 2)
                    < (winding_factor * FEM_data_matrix[i, param["window_h"]] * mc1.data_matrix[i, param["window_w"]])):
                continue

            # MagneticComponent class object
            geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, silent=True)

            core = fmt.Core(core_inner_diameter=FEM_data_matrix[i, param["core_w"]], window_w=FEM_data_matrix[i, param["window_w"]],
                            window_h=FEM_data_matrix[i, param["window_h"]],
                            # material="N95", temperature=25, frequency=freq, datasource="manufacturer_datasheet")
                            # material="95_100")
            mu_rel=3000, phi_mu_deg=10,
            sigma=0.5)
            # TODO: (completed)
            # mu_rel=3000, phi_mu_deg=10,
            # sigma=0.5)
            geo.set_core(core)

            # 3. set air gap parameters
            air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
            if int(FEM_data_matrix[i, param["n_air_gaps"]]) == 1:
                air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_h"]],
                                     FEM_data_matrix[i, param["air_gap_position"]])
            else:
                if int(FEM_data_matrix[i, param["mult_air_gap_type"]]) == 1:
                    position_list = list(np.linspace(0, 100, int(FEM_data_matrix[i, param["n_air_gaps"]])))
                    for position in position_list:
                        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_h"]],
                                             position)

                elif int(FEM_data_matrix[i, param["mult_air_gap_type"]]) == 2:
                    position_list = list(np.linspace(0, 100, int(FEM_data_matrix[i, param["n_air_gaps"]]) + 2))
                    position_list.remove(0.0)
                    position_list.remove(100.0)
                    for position in position_list:
                        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, FEM_data_matrix[i, param["air_gap_h"]],
                                             position)
            geo.set_air_gaps(air_gaps)

            # 4. set insulations
            insulation = fmt.Insulation()
            insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
            insulation.add_winding_insulations([0.0005], 0.0001)
            geo.set_insulation(insulation)

            # 5. create winding window and virtual winding windows (vww)
            winding_window = fmt.WindingWindow(core, insulation)
            vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

            # 6. create conductor and set parameters: use solid wires
            winding = fmt.Conductor(0, fmt.Conductivity.Copper)
            if j < len(litz_conductor_r):
                winding.set_litz_round_conductor(conductor_radius=litz_conductor_r[j], number_strands=litz_strand_n[j],
                                                 strand_radius=litz_strand_r[j], fill_factor=None,
                                                 conductor_arrangement=fmt.ConductorArrangement.Square)
            else:
                winding.set_solid_round_conductor(conductor_radius=conductor_r_list[j],
                                                  conductor_arrangement=fmt.ConductorArrangement.Square)
            # winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,
            # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

            # 7. add conductor to vww and add winding window to MagneticComponent
            vww.set_winding(winding, int(FEM_data_matrix[i, param["no_of_turns"]]), None)
            geo.set_winding_window(winding_window)

            try:
                # 5. create the model
                geo.create_model(freq=freq, visualize_before=False, save_png=False)

                # 6. start simulation
                geo.single_simulation(freq=freq, current=[i_max], show_results=False)

                shutil.copy2(src_path, os.path.join(os.path.dirname(__file__), save_directory_name))
                old_filename = os.path.join(os.path.dirname(__file__), save_directory_name, "log_electro_magnetic.json")
                new_filename = os.path.join(os.path.dirname(__file__), save_directory_name, f"case{counter3}.json")
                os.rename(old_filename, new_filename)
                working_directories.append(new_filename)
                file_names.append(f"case{counter3}")
                counter3 = counter3 + 1
                print(f"{counter3} of {n_cases_FEM * 2}")

            except (Exception,) as e:
                print("next iteration")
                logging.exception(e)

            # # 5. create the model
            # geo.create_model(freq=freq, visualize_before=False, save_png=False)
            #
            # # 6. start simulation
            # geo.single_simulation(freq=freq, current=[i_max], show_results=False)


def load_design(load_directory_name):
    working_directories = []
    labels = []
    working_directory = os.path.join(os.path.dirname(__file__), load_directory_name)
    print("##########################")
    print(f"{working_directory =}")
    print("##########################")
    file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]
    file_names.sort()
    counter2 = 0
    for name in file_names:
        temp_var = os.path.join(os.path.dirname(__file__), load_directory_name, name)
        working_directories.append(temp_var)
        # labels.append(f"case{counter2}")
        labels.append(name)
        counter2 = counter2 + 1

    zip_iterator = zip(file_names, working_directories)
    logs = dict(zip_iterator)

    # After the simulations the sweep can be analyzed
    # This could be done using the FEMMTLogParser:
    log_parser = fmt.FEMMTLogParser(logs)

    # In this case the self inductivity of winding1 will be analyzed
    inductivities = []
    active_power = []
    total_volume = []
    total_cost = []
    for name, data in log_parser.data.items():
        inductivities.append(data.sweeps[0].windings[0].self_inductance)
        active_power.append(data.sweeps[0].windings[0].active_power)
        total_volume.append(data.core_2daxi_total_volume)
        total_cost.append(data.total_cost)

    real_inductance = []
    for i in range(len(active_power)):
        real_inductance.append(inductivities[i].real)
    print(file_names)
    print(total_volume)
    print(active_power)
    print(total_cost)

    names = np.array(labels)
    c = np.random.randint(1, 5, size=len(active_power))

    norm = plt.Normalize(1, 4)
    cmap = plt.cm.RdYlGn

    fig, ax = plt.subplots()
    fmt.plt.title("Volume vs total losses")
    fmt.plt.xlabel("Volume (in cubic m)")
    fmt.plt.ylabel("Total losses (in W)")
    sc = plt.scatter(total_volume, active_power, c=c, s=50, cmap=cmap, norm=norm)

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        # text = "{}, {}".format(" ".join(list(map(str, ind["ind"]))),
        #                       " ".join([names[n] for n in ind["ind"]]))
        text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    ax.grid()
    plt.show()


def compare_graph(directory_names):
    temp = 0

    for directory_name in directory_names:
        working_directories = []
        labels = []
        working_directory = os.path.join(os.path.dirname(__file__), directory_name)
        print("##########################")
        print(f"{working_directory =}")
        print("##########################")
        file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]
        file_names.sort()
        counter2 = 0
        for name in file_names:
            temp_var = os.path.join(os.path.dirname(__file__), directory_name, name)
            working_directories.append(temp_var)
            # labels.append(f"case{counter2}")
            labels.append(name)
            counter2 = counter2 + 1

        zip_iterator = zip(file_names, working_directories)
        logs = dict(zip_iterator)

        # After the simulations the sweep can be analyzed
        # This could be done using the FEMMTLogParser:
        log_parser = fmt.FEMMTLogParser(logs)

        # In this case the self inductivity of winding1 will be analyzed
        inductivities = []
        active_power = []
        total_volume = []
        for name, data in log_parser.data.items():
            inductivities.append(data.sweeps[0].windings[0].self_inductance)
            active_power.append(data.sweeps[0].windings[0].active_power)
            total_volume.append(data.core_2daxi_total_volume)

        real_inductance = []
        for i in range(len(active_power)):
            real_inductance.append(inductivities[i].real)

        print(real_inductance)
        print(active_power)

        # fig, ax = fmt.plt.subplots()  # Create a figure containing a single axes.
        # fmt.plt.title("Normalised volume vs Normalised losses")
        # fmt.plt.xlabel("Normalised Volume")
        # fmt.plt.ylabel("Normalised Losses")
        # ax.grid()
        if temp == 0:
            plt.plot(total_volume, active_power, 'o')
            temp =1
        else:
            plt.plot(total_volume, active_power, '*')
    plt.title("volume vs losses")
    plt.xlabel("Total Volume")
    plt.ylabel("Total Losses")
    plt.grid()
    fmt.plt.show()


def plot_3d(load_directory_name):
    working_directories = []
    labels = []
    working_directory = os.path.join(os.path.dirname(__file__), load_directory_name)
    print("##########################")
    print(f"{working_directory =}")
    print("##########################")
    file_names = [f for f in listdir(working_directory) if isfile(join(working_directory, f))]
    print(file_names)
    num = [re.findall(r'\d+', item) for item in file_names]
    case_num_list = [int(item[0]) for item in num]
    file_names.sort()
    print(file_names)
    counter2 = 0
    for name in file_names:
        temp_var = os.path.join(os.path.dirname(__file__), load_directory_name, name)
        working_directories.append(temp_var)
        # labels.append(f"case{counter2}")
        labels.append(name)
        counter2 = counter2 + 1

    zip_iterator = zip(file_names, working_directories)
    logs = dict(zip_iterator)

    # After the simulations the sweep can be analyzed
    # This could be done using the FEMMTLogParser:
    log_parser = fmt.FEMMTLogParser(logs)

    # In this case the self inductivity of winding1 will be analyzed
    inductivities = []
    active_power = []
    total_volume = []
    total_cost = []
    for name, data in log_parser.data.items():
        inductivities.append(data.sweeps[0].windings[0].self_inductance)
        active_power.append(data.sweeps[0].windings[0].active_power)
        total_volume.append(data.core_2daxi_total_volume)
        total_cost.append(data.total_cost)

    real_inductance = []
    for i in range(len(active_power)):
        real_inductance.append(inductivities[i].real)
    print(file_names)
    print(total_volume)
    print(active_power)
    print(total_cost)

    names = np.array(labels)
    case_num_list = np.array(case_num_list)

    X = np.zeros((len(total_volume), 1))
    X = np.hstack((X, np.reshape(total_volume, (len(total_volume), 1))))
    X = np.hstack((X, np.reshape(active_power, (len(active_power), 1))))
    X = np.hstack((X, np.reshape(total_cost, (len(total_cost), 1))))
    X = np.hstack((X, np.reshape(case_num_list, (len(case_num_list), 1))))
    X = np.delete(X, 0, 1)
    X = X[X[:, 3].argsort()]
    print(len(X))
    return X, names


def visualize3DData (K):
    names = K[1]
    X = K[0]
    """Visualize data in 3d plot with popover next to mouse position.

    Args:
        X (np.array) - array of points, of shape (numPoints, 3)
    Returns:
        None
    """
    fig = plt.figure(figsize = (16,10))
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_title('Volume vs loss vs cost')
    ax.set_xlabel('volume (in cubic m)')
    ax.set_ylabel('losses (in W)')
    ax.set_zlabel('cost (in Euro)')

    c = np.random.randint(1, 5, size=len(X))

    norm = plt.Normalize(1, 4)
    cmap = plt.cm.RdYlGn

    ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=c, s=10, cmap=cmap, norm=norm, depthshade = False, picker = True)

    def distance(point, event):
        """Return distance between mouse position and given data point

        Args:
            point (np.array): np.array of shape (3,), with x,y,z in data coords
            event (MouseEvent): mouse event (which contains mouse position in .x and .xdata)
        Returns:
            distance (np.float64): distance (in screen coords) between mouse pos and data point
        """
        assert point.shape == (3,), "distance: point.shape is wrong: %s, must be (3,)" % point.shape

        # Project 3d data space to 2d data space
        x2, y2, _ = proj3d.proj_transform(point[0], point[1], point[2], plt.gca().get_proj())
        # Convert 2d data space to 2d screen space
        x3, y3 = ax.transData.transform((x2, y2))

        return np.sqrt ((x3 - event.x)**2 + (y3 - event.y)**2)

    def calcClosestDatapoint(X, event):
        """"Calculate which data point is closest to the mouse position.

        Args:
            X (np.array) - array of points, of shape (numPoints, 3)
            event (MouseEvent) - mouse event (containing mouse position)
        Returns:
            smallestIndex (int) - the index (into the array of points X) of the element closest to the mouse position
        """
        distances = [distance (X[i, 0:3], event) for i in range(X.shape[0])]
        return np.argmin(distances)

    def annotatePlot(X, index):
        """Create popover label in 3d chart

        Args:
            X (np.array) - array of points, of shape (numPoints, 3)
            index (int) - index (into points array X) of item which should be printed
        Returns:
            None
        """
        # If we have previously displayed another label, remove it first
        if hasattr(annotatePlot, 'label'):
            annotatePlot.label.remove()
        # Get data point from array of points X, at position index
        x2, y2, _ = proj3d.proj_transform(X[index, 0], X[index, 1], X[index, 2], ax.get_proj())
        annotatePlot.label = plt.annotate( "case %d" % index,
            xy = (x2, y2), xytext = (-20, 20), textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        fig.canvas.draw()

    def onMouseMotion(event):
        """Event that is triggered when mouse is moved. Shows text annotation over data point closest to mouse."""
        closestIndex = calcClosestDatapoint(X, event)
        annotatePlot (X, closestIndex)

    fig.canvas.mpl_connect('motion_notify_event', onMouseMotion)  # on mouse motion
    plt.show()


if __name__ == '__main__':

    # automated_design_func()

    # design_name = "sweep1_31_10_22"
    # load_design(design_name)

    design_name = "sweep_new_3d_with_cost"
    K = plot_3d(design_name)
    print(K)
    visualize3DData(K)

    # directory_names = ["sweep_examples_5", "hyperbolic_filtered_0_01_0_6"]
    # compare_graph(directory_names)






























# def input_pre_process(manual_core_w: list, manual_window_h: list, manual_window_w: list, db_core_names: list, db_litz_names: list,
#                       material_names: list, mult_air_gap_type: list):
#     all_manual_combinations = list(product(manual_core_w, manual_window_h, manual_window_w))
#     manual_core_w = [item[0] for item in all_manual_combinations]
#     manual_window_h = [item[1] for item in all_manual_combinations]
#     manual_window_w = [item[2] for item in all_manual_combinations]
#
#     core_db = fmt.core_database()
#     db_core_w = [core_db[core_name]["core_inner_diameter"] for core_name in db_core_names]
#     db_window_h = [core_db[core_name]["window_h"] for core_name in db_core_names]
#     db_window_w = [core_db[core_name]["window_w"] for core_name in db_core_names]
#
#     core_w_list = db_core_w + manual_core_w
#     window_h_list = db_window_h + manual_window_h
#     window_w_list = db_window_w + manual_window_w
#
#     litz_db = fmt.litz_database()
#     litz_conductor_r = [litz_db[litz_name]["conductor_radii"] for litz_name in db_litz_names]
#     litz_strand_r = [litz_db[litz_name]["strand_radii"] for litz_name in db_litz_names]
#     litz_strand_n = [litz_db[litz_name]["strands_numbers"] for litz_name in db_litz_names]
#
#     # material_data_list = [material_db.permeability_data_to_pro_file(25, freq, material_name, "manufacturer_datasheet")
#     #                       for material_name in material_names]
#     #
#     # mu_imag_interpol_func = [interp1d(material_data_list[i][0], material_data_list[i][1], kind="cubic") for i in
#     #                          range(len(material_names))]
#
#
#     mu_rel = [material_db.get_material_property(material_name=material_name, property="initial_permeability")
#               for material_name in material_names]
#
#     mult_air_gap_type_list = []
#     for item in mult_air_gap_type:
#         if item == 'edge_distributed':
#             mult_air_gap_type_list.append(1)
#         elif item == 'center_distributed':
#             mult_air_gap_type_list.append(2)
#         else:
#             raise Exception('Wrong string input for multiple air-gap type')
#
#     return core_w_list, window_h_list, window_w_list, litz_conductor_r, litz_strand_r, litz_strand_n, mu_rel, mult_air_gap_type_list
#
#
# def filter_1(data_matrix, param, litz_conductor_r, solid_conductor_r):
#     min_conductor_r = min(litz_conductor_r + solid_conductor_r)
#     data_matrix = data_matrix[np.where((data_matrix[:, param["no_of_turns"]] * np.pi * min_conductor_r ** 2)
#                                                  < (winding_factor * data_matrix[:,
#                                                                      param["window_h"]] * data_matrix[:,
#                                                                                           param[
#                                                                                               "window_w"]]))]
#     return data_matrix
#
#
# def filter_2(data_matrix, param, inductance_percent_tolerance, goal_inductance):
#
#     data_matrix = data_matrix[np.where((data_matrix[:, param["inductance"]] > ((100 - inductance_percent_tolerance) / 100)
#                                             * goal_inductance) & (data_matrix[:, param["inductance"]] <
#                                             ((100 + inductance_percent_tolerance) / 100) * goal_inductance))]
#     return data_matrix
#
#
# def filter_3(data_matrix, param, material_names, percent_of_B_sat):
#     B_sat_dict = {}
#     counter1 = 0
#     for material_name in material_names:
#         B_sat_key = material_db.get_material_property(material_name=material_name, property="initial_permeability")
#         B_sat_dict[B_sat_key] = material_db.get_material_property(material_name=material_name, property="max_flux_density")
#         counter1 = counter1 + 1
#
#     # Creating B_saturated array corresponding to the material type
#     B_sat = np.zeros((len(data_matrix), 1))
#     for index in range(len(data_matrix)):
#         B_sat[index] = B_sat_dict[data_matrix[index, param["mu_rel"]]]
#
#     # flux_max = L * i_max / N
#     total_flux_max = (data_matrix[:, param["inductance"]] * current_max) / data_matrix[:, param["no_of_turns"]]
#     data_matrix_2 = np.hstack((data_matrix, np.reshape(total_flux_max, (len(total_flux_max), 1))))
#     param["total_flux_max"] = 15
#
#     B_max_center = total_flux_max / data_matrix_2[:, param["center_leg_area"]]
#     B_max_middle = total_flux_max / (
#             np.pi * data_matrix_2[:, param["core_w"]] * data_matrix_2[:, param["core_h_middle"]])
#     B_max_outer = total_flux_max / data_matrix_2[:, param["outer_leg_area"]]
#
#     data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_center, (len(B_max_center), 1))))
#     param["B_max_center"] = 16
#     data_matrix_2 = np.hstack((data_matrix_2, np.reshape(B_max_outer, (len(B_max_outer), 1))))
#     param["B_max_outer"] = 17
#
#     data_matrix_3 = np.zeros((0, 18))
#     for index in range(len(data_matrix_2)):
#         if (data_matrix_2[index, param["B_max_center"]] < (percent_of_B_sat / 100) * B_sat[index]) & \
#                 (data_matrix_2[index, param["B_max_outer"]] < (percent_of_B_sat / 100) * B_sat[index]):
#             data_matrix_3 = np.vstack([data_matrix_3, data_matrix_2[index, :]])
#     return data_matrix
#
#
# # User inputs (DESIGN PARAMETERS)--------------------------------------------------------------------------------------
# save_directory_name = "sweep_new_3d_with_cost"
# component = 'inductor'
# goal_inductance = 120 * 1e-6
# inductance_percent_tolerance = 10       # Inductance tolerance of +-10% is applied
# winding_factor = 0.91
# current_max = 8                         # Max current amplitude with assumption of sinusoidal current waveform
# percent_of_b_sat = 70                   # Percent of B_sat allowed in the designed core
# percent_of_total_loss = 30              # Percent of total_loss allowed in FEM simulation
# freq = 100 * 1e3                        # Switching frequency in Hz
# copper_conductivity = 5.96 * 1e7        # copper conductivity (sigma) @ 20 degree celsius (in siemens/meter)
#
# # Set core-geometry from core database or/and manual entry
# db_core_names = []  # "PQ 40/40", "PQ 40/30"
#
# manual_core_w = list(np.linspace(0.005, 0.05, 10))
# manual_window_h = list(np.linspace(0.01, 0.08, 5))
# manual_window_w = list(np.linspace(0.005, 0.04, 10))
#
# # Set winding settings (Solid and Litz winding type)
# db_litz_names = ["1.5x105x0.1"]  # "1.5x105x0.1", "1.4x200x0.071"
# solid_conductor_r = [0.0013]
#
# # Set air-gap and core parameters
# no_of_turns = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]  # Set No. of turns (N)
# n_air_gaps = [1, 2]                                     # Set No. of air-gaps (n)
# air_gap_height = list(np.linspace(0.0001, 0.0005, 5))   # Set air-gap height in metre (l)
# air_gap_position = list(np.linspace(20, 80, 2))         # Set air-gap position in percent w.r.t. core window height
#
# # Set core material
# material_list = ["N95"]
#
# # Multiple air-gap type ('edge_distributed' and/or 'center_distributed')
# mult_air_gap_type = ['center_distributed']
#
# # ---------------------------------------------------------------------------------------------------------------------
# core_w_list, window_h_list, window_w_list, litz_conductor_r, litz_strand_r, litz_strand_n, mu_rel, mult_air_gap_type_list = \
#     input_pre_process(manual_core_w, manual_window_h, manual_window_w, db_core_names, db_litz_names, material_list, mult_air_gap_type)
#
# # Call to Reluctance model (Class MagneticCircuit)
# mc1 = fmt.MagneticCircuit(core_w=core_w_list, window_h=window_h_list, window_w=window_w_list, no_of_turns=no_of_turns,
#                           n_air_gaps=n_air_gaps, air_gap_h=air_gap_height, air_gap_position=air_gap_position,
#                           mu_rel=mu_rel, mult_air_gap_type=mult_air_gap_type_list, air_gap_method='percent',
#                           component_type=component, sim_type='sweep')
# param = mc1.get_param_pos_dict()
#
# # Filtration of the design cases which are not important
# data_matrix_1 = filter_1(mc1.data_matrix, param, litz_conductor_r, solid_conductor_r)
# data_matrix_2 = filter_2(data_matrix_1, param, inductance_percent_tolerance, goal_inductance)
# data_matrix_3 = filter_3(data_matrix_2, param, litz_conductor_r, solid_conductor_r)
# #
# # data_matrix_4 = filter_4(data_matrix_3, param, litz_conductor_r, solid_conductor_r)
# #
# # fem_simulation()