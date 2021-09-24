import FEMMT
from FEMMT import MagneticComponent, r_basis, sigma, r_cyl_cyl, r_cheap_cyl_cyl, r_round_inf
import numpy as np
from scipy.optimize import fsolve, brentq


# TODO: sicheres & schnelles Ausschließen der geometrisch unmöglichen Konfigurationen auf high level Ebene
# TODO: bis Ende der Woche ein Streupfad Reluktanzmodell "in Reihe" mit einer Induktivitätssimulation und
#  Verlustsimulation laufen haben


#                                               -- Definitions --
# ----------------------------------------------------------------------------------------------------------------------
mu0 = 4e-7*np.pi
mur = 3000
inductance_calculation = 1
skin_accuracy = 0.5

# Goal Parameters
L_s1 = 6.8e-5
L_h = 10 * L_s1
n = 3.2

# Fixed Parameters
frequency = 250000
window_h = 0.0295
window_w = 0.012
core_w = 0.015

# Heuristic air gap lengths
# air_gap_h_bot = 0.0005
# air_gap_h_top = 0.0005
# air_gap_h_stray = 0.0005

# Sweep Parameters
N1 = 27  # Turns in main window
N2 = 7  # Turns in main window
Ns1 = 5  # Turns in stray window
Ns2 = 6  # Turns in stray window
lower_air_gap_pos = 30
strand_radius = 0.025e-3
N_strands_prim = 300
N_strands_sec = 450


#iter..
#geo = MagneticComponent
#geo.relu_pre(iter.., goals, core_paras)





#                                           -- Reluctance Model --
# ----------------------------------------------------------------------------------------------------------------------

# Inductance Matrix
L_11 = L_s1 + L_h
L_22 = L_h / n**2
M = L_h / n
L = np.array([[L_11, M], [M, L_22]])

# Winding Matrix
N = np.array([[N1, N2], [Ns1, Ns2]])

# Reluctance Matrix
R = np.matmul(np.matmul(N, np.linalg.inv(L)), np.transpose(N))

# Goal Reluctance Values
R_stray = -R[0, 1]
R_top = R[0, 0] - R_stray
R_bot = R[1, 1] - R_stray

# Print Reluctances
print(f"\n-------------------------------\n"
      f"Goal Values:\n"
      f"R = {R}\n"
      f"R_stray = {R_stray}\n"
      f"R_top = {R_top}\n"
      f"R_bot = {R_bot}\n")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Idealized Air gaps [needed for initializing the Newton Solver]
# - Reluctance of air gap
air_gap_h_stray_0 = R_stray * mu0 * 2 * (core_w / 2 + window_w) * np.pi * (window_h * (50 - lower_air_gap_pos) / 100)
# - Reluctance of core and air gap
A_core = (core_w / 2) ** 2 * np.pi
l_core_top = 2 * 0.5 * window_h + window_w
air_gap_h_top_0 = (R_top - l_core_top / (mu0 * mur * A_core)) * (mu0 * (A_core))
# - Reluctance of core and air gap
l_core_bot = 2 * (lower_air_gap_pos / 100) * window_h + window_w
air_gap_h_bot_0 = (R_bot - l_core_bot / (mu0 * mur * A_core)) * mu0 * A_core


# Reluctances calculated based on [A Novel Approach for 3D Air Gap Reluctance Calculations]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def R_top_a(air_gap_h_top):
    return r_round_inf(l=air_gap_h_top, r=core_w/2,
                       sigma=sigma(l=air_gap_h_top, w=core_w/2, R_equivalent=r_basis(l=air_gap_h_top,
                                                                                     w=core_w, h=window_h*0.5))
                       ) - R_top


air_gap_h_top = brentq(R_top_a, 1e-6, window_h/4)
print(f"With fringing: air_gap_h_top = {air_gap_h_top}\n"
      f"Idealized: air_gap_h_top_0 = {air_gap_h_top_0}\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def R_bot_a(air_gap_h_bot):
    return r_round_inf(l=air_gap_h_bot, r=core_w/2,
                       sigma=sigma(l=air_gap_h_bot, w=core_w/2, R_equivalent=r_basis(l=air_gap_h_bot, w=core_w,
                                                                                     h=window_h*lower_air_gap_pos/100))
                       ) - R_bot


air_gap_h_bot = brentq(R_bot_a, 1e-6, window_h/4)
print(f"With fringing: air_gap_h_bot = {air_gap_h_bot}\n"
      f"Idealized: air_gap_h_bot_0 = {air_gap_h_bot_0}\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#R_stray_cheap = R_cheap_cyl_cyl(r_o=window_w+core_w/2, l=air_gap_h_stray, w=w_stray)
def R_stray_a(air_gap_h_stray):
    w_stray = window_h * (50 - lower_air_gap_pos) / 100 - air_gap_h_bot / 2 - air_gap_h_top_0 / 2
    return r_cyl_cyl(l=air_gap_h_stray,
                    sigma=sigma(l=air_gap_h_stray,
                                w=w_stray,
                                R_equivalent=0.5*r_basis(l=air_gap_h_stray,
                                                         w=w_stray,
                                                         h=window_w-air_gap_h_stray)),
                    w=w_stray,
                    r_o=window_w+core_w/2) - R_stray


air_gap_h_stray = brentq(R_stray_a, 1e-6, window_w/2)
print(f"With fringing: air_gap_h_stray = {air_gap_h_stray}\n"
      f"Idealized: air_gap_h_stray_0 = {air_gap_h_stray_0}\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# R_a = np.array([[R_stray+R_top, R_stray], [R_stray, R_stray+R_bot]])  # print(f"\nN = {N}")
# L_a = np.matmul(np.matmul(np.transpose(N), np.linalg.inv(R_a)), N)
# print(f"\nN = {N}\n"
#       f"L_a= {L_a}\n")





# CHECK Air Gap Data
if (air_gap_h_stray<=0) or (air_gap_h_top<=0) or (air_gap_h_bot<=0):
    print(f"Negative Air gap would be needed to realize desired Inductance value!")
else:
    #                                           -- CHECK Saturation --
    # ------------------------------------------------------------------------------------------------------------------
    do_simulate = True

    if do_simulate:      # TODO: Sättigung überprüfen

        #                                         -- FEM Simulation --
        # --------------------------------------------------------------------------------------------------------------

        # -- FFT --
        I0 = 6 * 1.25 # 6 is peak of current wave_form
        # -- Component Design --
        geo = MagneticComponent(component_type="transformer")

        geo.update_core(core_type="EI", window_h=window_h, window_w=window_w, core_w=core_w,
                        non_linear=False, core_material = 95_100, core_re_mu_rel=mur)

        geo.update_air_gaps(method="percent",
                            n_air_gaps=2,
                            position_tag=[0, 0],
                            air_gap_h=[air_gap_h_bot, air_gap_h_top],
                            air_gap_position=[lower_air_gap_pos, 50],
                            stray_path=[[0, 1], geo.core_w/2+geo.window_w-air_gap_h_stray])

        geo.update_conductors(n_turns=[[N1, Ns1], [N2, Ns2]],
                              conductor_type=["litz", "litz"],
                              winding=["interleaved", "interleaved"],
                              scheme=["horizontal", "horizontal"],  # ["square", "square"]  "horizontal"
                              conductor_radii=[0.0005, 0.0009],
                              litz_para_type=['implicite_litz_radius', 'implicite_litz_radius'],
                              strands_numbers=[N_strands_prim, N_strands_sec],
                              ff=[0.5, 0.5],
                              strand_radii=[strand_radius, strand_radius],
                              cond_cond_isolation=[0.0002, 0.0002, 0.0004],
                              core_cond_isolation=[0.002])

        # -- Simulation --
        if inductance_calculation:
            geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=skin_accuracy)

            print(f"\n-------------------------------\n"
                  f"Reluctance Model Values:\n"
                  f"L_11 = {L_11}\n"
                  f"L_22 = {L_22}\n"
                  f"M = {M}\n"

                  f"R = {R}\n"
                  f"R_stray_a = {R_stray}\n"
                  f"R_top_a = {R_top}\n"
                  f"R_bot_a = {R_bot}\n"
                  # f"R_stray_cheap = {R_stray_cheap}\n"
                  f"\n"
                  # f"air_gap_h_stray = {air_gap_h_stray}\n"
                  # f"air_gap_h_top = {air_gap_h_top}\n"
                  # f"air_gap_h_bot = {air_gap_h_bot}\n"
                  )

            L_FEM_11 = geo.L_h_conc + geo.L_s_conc  # TODO:Read Mühlethaler Diss Bereich sinnvoll approx w/h
            L_FEM_22 = geo.L_h_conc / geo.ü_conc ** 2
            M_FEM = geo.L_h_conc / geo.ü_conc

            L_FEM = np.array([[L_FEM_11, -M_FEM], [-M_FEM, L_FEM_22]])
            """
            print(f"\n"
                  f"L_FEM = {L_FEM}\n"
                  f"L_FEM-1 = {np.linalg.inv(L_FEM)}\n"
                  f"det(L_FEM-1) = {np.linalg.det(np.linalg.inv(L_FEM))}\n"
                  f"L_a = {L_a}\n"
                  f"L_a-1 = {np.linalg.inv(L_a)}\n"
                  f"det(L_a-1) = {np.linalg.det(np.linalg.inv(L_a))}\n"
                  f"L_a-1/4 = {np.linalg.inv(L_a)/4}")
            """

            R_FEM = np.matmul(np.matmul(N, np.linalg.inv(L_FEM)), np.transpose(N))
            print(f"TEST:{np.matmul(np.linalg.inv(L_FEM), L_FEM)}")
            # R_something = np.matmul(np.matmul(N, np.linalg.inv(L_a)), np.transpose(N))
            # print(L_FEM)
            R_FEM_stray = -R_FEM[0, 1]
            R_FEM_top = R_FEM[0, 0] - R_FEM_stray
            R_FEM_bot = R_FEM[1, 1] - R_FEM_stray
            print(f"\n-------------------------------\n"
                  f"FEM Values:\n"
                  f"L_FEM_11 = {L_FEM_11}\n"
                  f"L_FEM_22 = {L_FEM_22}\n"
                  f"M_FEM = {M_FEM}\n"
                  # f"R_something = {R_something}\n"
                  f"R_FEM = {R_FEM}\n"
                  f"R_FEM_stray = {R_FEM_stray}\n"
                  f"R_FEM_top = {R_FEM_top}\n"
                  f"R_FEM_bot = {R_FEM_bot}\n")

        else:
            geo.single_simulation(freq=frequency, current=[I0, N1/N2*I0], phi=[0, 180], skin_mesh_factor=skin_accuracy)


    else:
        print(f"Core is saturated!")


"""
# -- Sweep Parameters --
# - Windings
N1 = 27 # Turns in main window
N2 = 7  # Turns in main window
Ns1 = 5  # Turns in stray window
Ns2 = 6  # Turns in stray window
# - Air gaps
air_gap_h_bot = 0.0012
air_gap_h_top = 0.0002
air_gap_h_stray = 0.0002
# - Core
lower_air_gap_pos = 30
# - Conductors
strand_radius = 0.025e-3
N_strands_prim = 400
N_strands_sec = 800
"""