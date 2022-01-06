from femmt import MagneticComponent

# Create Object
geo = MagneticComponent(component_type="integrated_transformer")

# Goal Parameters
L_s1 = 6.8e-5
L_h = 10 * L_s1
n = 3.2

# Fixed Parameters
frequency = 250000
window_h = 0.0295
window_w = 0.012
core_w = 0.015

# Sweep Parameters
N1 = 27  # Turns in main window
N2 = 7  # Turns in main window
Ns1 = 5  # Turns in stray window
Ns2 = 6  # Turns in stray window
midpoint = 30
strand_radius = 0.025e-3
N_strands_prim = 300
N_strands_sec = 450

# Inductance Matrix
L_11 = L_s1 + L_h
L_22 = L_h / n**2
M = L_h / n
L = [[L_11, M],
     [M, L_22]]

# 720941 oben
# 3179294 unten
# 328117 mitte



# Winding Matrix
N1 = [[N1, N2],
     [Ns1, Ns2]]

# List of lists
N_init = [N1]
Core_init = [{'window_h': window_h, 'window_w': window_w, 'core_w': core_w}]
Stray_init = [{'midpoint': (midpoint - 50) / 100 * window_h}]

#geo.core.update(type="EI", window_h=window_h, window_w=window_w, core_w=core_w,
#                           non_linear=False, material=95_100, re_mu_rel=3000)

res = geo.reluctance_model.air_gap_design(L_goal=L,
                                          N_init=N_init, core_par_init=Core_init, stray_path_init=Stray_init,
                                          current=[7.5, -7.5*3.2], Bmax=0.05)

print(res)


# einfaches LTspice Modell zum Übersetzungsverhältnis