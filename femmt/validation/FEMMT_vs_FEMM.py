import femmt as fmt
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import leapythontoolbox as lpt
import os

if not os.path.isdir("/results"):
    os.mkdir("/results")

def inductor_1():
    geo = fmt.MagneticComponent(component_type="inductor")

    frequency = 250000

    geo.global_accuracy = 0.3

    # Update Geometry
    geo.core.update(window_h=0.03, window_w=0.011)

    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.002], position_tag=[0])

    geo.update_conductors(n_turns=[[26]], conductor_type=["solid"], conductor_radii=[0.001],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0001],
                          conductivity_sigma=["copper"])

    geo.create_model(freq=frequency, visualize_before=False, do_meshing=True, save_png=True)

    geo.single_simulation(freq=frequency, current=[2.5], show_results=True)

def inductor_2():
    geo = fmt.MagneticComponent(component_type="inductor")

    geo.global_accuracy = 0.3

    frequency = 250000
    # Update Geometry
    geo.core.update(window_h=0.03, window_w=0.011)

    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.001], position_tag=[0])

    geo.update_conductors(n_turns=[[13]], conductor_type=["solid"], conductor_radii=[0.001],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0001],
                          conductivity_sigma=["copper"])

    geo.create_model(freq=frequency, visualize_before=False, do_meshing=True, save_png=True)

    geo.single_simulation(freq=frequency, current=[2.5], show_results=False)

def inductor_3():
    geo = fmt.MagneticComponent(component_type="inductor")
    frequency = 250000
    geo.global_accuracy = 0.1

    # Update Geometry
    geo.core.update(window_h=0.03, window_w=0.011)

    geo.air_gaps.update(method="center", n_air_gaps=1, air_gap_h=[0.001], position_tag=[0])

    geo.update_conductors(n_turns=[[20]], conductor_type=["solid"], conductor_radii=[0.001],
                          winding=["primary"], scheme=["square"],
                          core_cond_isolation=[0.001, 0.001, 0.002, 0.001], cond_cond_isolation=[0.0001],
                          conductivity_sigma=["copper"])

    geo.create_model(freq=frequency, visualize_before=False, do_meshing=True, save_png=True)

    geo.single_simulation(freq=frequency, current=[2.5], show_results=False)

def femmt_ind_loss_current_vector(inductor_vector: str, current_vector):

    loss_vector = np.zeros([len(inductor_vector), len(current_vector)])
    current_vector_str = []
    for count, current in enumerate(current_vector):
        string = str(current)
        current_string = ""
        for letter in string:
            if letter == ".":
                current_string += "a"
            else:
                current_string += letter

        current_vector_str.append(current_string)
    print(current_vector_str)

    for p in Path('.').glob('*.json'):
        for count_inductor, inductor in enumerate(inductor_vector):
            for count_current, current in enumerate(current_vector_str):
                f = open(f'femmt_{inductor}_{current}.json')
                data = json.load(f)
                loss_vector[count_inductor][count_current] = data["Losses"]["all_windings"]


    return loss_vector
        #print(f"{p.name}:\n{p.read_text()}\n")

def femmt_ind_loss_frequency_vector(inductor: str, current: str, frequency_list: list):
    loss_vector = []

#    for p in Path('.').glob('*.json'):

    for frequency in frequency_list:
        f = open(f'femmt_{inductor}_{frequency}kHz_{current}.json')
        data = json.load(f)
        #print(data["Losses"]["all_windings"])
        loss_vector.append(data["Losses"]["all_windings"])

    return loss_vector

def plot_results():

    frequency_sweep = True

    if frequency_sweep:
        plt.figure(figsize=(4, 2.5))
        frequency_vector = [0, 50, 100, 150, 200, 250]

        # FEMM
        femm_ind1_frequency_loss_vector = [0.0786535, 3.39698,  5.26505, 6.76725,  8.06984, 9.24349] #, 10.3224]
        femm_ind2_frequency_loss_vector = [0.0363803,  0.746134,  1.212,  1.59707,  1.93515,  2.2426] #,  2.52729]
        femm_ind3_frequency_loss_vector = [0.0591427,  2.17183, 3.40977,  4.41266,  5.2855,  6.07401]

        # FEMMT
        femmt_ind1_frequency_loss_vector = femmt_ind_loss_frequency_vector("ind1", "2a5", frequency_vector)
        femmt_ind2_frequency_loss_vector = femmt_ind_loss_frequency_vector("ind2", "2a5", frequency_vector)
        femmt_ind3_frequency_loss_vector = femmt_ind_loss_frequency_vector("ind3", "2a5", frequency_vector)

        print(femmt_ind1_frequency_loss_vector)

        plt.plot(frequency_vector, femmt_ind1_frequency_loss_vector, 'o-', label="FEMMT", linewidth=2,
                 color=lpt.gnome_colors["blue"])
        plt.plot(frequency_vector, femm_ind1_frequency_loss_vector, 'x--', linewidth=2, color=lpt.gnome_colors["red"],
                 label="FEMM")
        plt.plot(frequency_vector, femmt_ind2_frequency_loss_vector, 'o-', linewidth=2,
                 color=lpt.gnome_colors["blue"])
        plt.plot(frequency_vector, femm_ind2_frequency_loss_vector, 'x--', linewidth=2, color=lpt.gnome_colors["red"])
        plt.plot(frequency_vector, femmt_ind3_frequency_loss_vector, 'o-', linewidth=2,
                 color=lpt.gnome_colors["blue"])
        plt.plot(frequency_vector, femm_ind3_frequency_loss_vector, 'x--', linewidth=2, color=lpt.gnome_colors["red"])

        x_koordinate = 260
        plt.text(x_koordinate, 9, r"$L_{\mathrm{1}}$")
        plt.text(x_koordinate, 2.3, r"$L_\mathrm{2}$")
        plt.text(x_koordinate, 6, r"$L_\mathrm{3}$")
        plt.xlim((0, 300))
        fs = 11
        plt.xlabel(r"$f$ / kHz", fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.ylabel(r"$P_{\mathrm{winding}}$ / W", fontsize=fs)
        plt.yticks(fontsize=fs)

        plt.grid()
        plt.legend()
        plt.xticks([0, 50, 100, 150, 200, 250, 300])
        plt.savefig("/results/losses_solid.pdf", bbox_inches="tight")
        plt.show()



    else:
        # Current sweep

        femm_100kHz_current_vector = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
        femm_ind1_loss_vector = [0.210602, 0.842408, 1.89542, 3.36963, 5.26505, 7.58167, 10.3195]

        femmt_ind1_loss_vector = femmt_ind_loss_current_vector(["ind1"], femm_100kHz_current_vector)
        print(np.shape(femm_100kHz_current_vector))
        print(np.shape(femmt_ind1_loss_vector))
        femmt_ind1_loss_vector.resize(7)

        plt.plot(femm_100kHz_current_vector, femm_ind1_loss_vector, label="femm")
        plt.plot(femm_100kHz_current_vector, femmt_ind1_loss_vector, label="femmt")
        plt.xlabel('Current in A')
        plt.ylabel('Losses in W')
        plt.legend()
        plt.grid()
        plt.show()




if __name__ == '__main__':
    #plot_results()
    inductor_2()
    inductor_1()
    inductor_3()

