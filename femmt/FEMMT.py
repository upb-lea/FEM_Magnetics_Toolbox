import numpy as np

# General properties
core_type = "EI"  # Basic shape of magnetic conductor
#conductor_type = "stacked"  # Vertical packing of conductors
#conductor_type = "full"  # One massive Conductor in each window
conductor_type = "foil"  # Horizontal packing of conductors
y_symmetric = 1  # Mirror-symmetry across y-axis
axi_symmetric = 1  # Axial-symmetric model (idealized full-cylindrical)
# Further geometry settings
s = 0.4  # Parameter for mesh-accuracy
n_air_gaps = 3  # Number of air gaps - Needs a re-visit for non axi-symmetric case
n_conductors = 6  # Number of (homogenised) conductors in one window
core_cond_isolation = 0.001  # gap between Core and Conductor
cond_cond_isolation = 0.0002  # gap between Core and Conductor
# Dimensions in meter
core_w = 0.02  # Axi symmetric case | core_w := core radius
window_w = 0.01  # Winding window width
window_h = 0.03  # Winding window height


# Pre-Settings
if core_type == "EI":
    n_windows = 2
# Symmetry: Choose between asymmetric, symmetric and axi symmetric
if y_symmetric == 0:
    axi_symmetric = 0
if axi_symmetric == 1:
    y_symmetric = 1
# Air Gap Data (random air gap generation)
if y_symmetric == 0:
    raise("Not up to date! Air Gap generation must be adopted from symmetric case")
    air_gaps = np.empty((n_air_gaps, 4))
    for i in range(0, n_air_gaps):
        if i < n_air_gaps/2:
            position_tag = -1
        else:
            position_tag = 1
        airgap_position = np.random.rand(1)*window_h-window_h/2
        airgap_h = np.random.rand(1)*0.005
        c_airgap = airgap_h / 3 * s
        air_gaps[i, :] = np.array([position_tag, airgap_position, airgap_h, c_airgap])

if y_symmetric == 1:
    air_gaps = np.empty((n_air_gaps, 4))
    i = 0
    while i in range(0, n_air_gaps):
        position_tag = 0
        airgap_h = 0.001#np.random.rand(1)*0.005 + 0.001
        airgap_position = np.random.rand(1)*(window_h-airgap_h)-(window_h/2-airgap_h/2)
        c_airgap = airgap_h / 3 * s
        # Overlapping Control
        for j in range(0, air_gaps.shape[0]):
            if position_tag == air_gaps[j, 0] and air_gaps[j, 1] + air_gaps[j, 2] / 2 > airgap_position > air_gaps[j, 1] - air_gaps[j, 2] / 2:
                print("Overlapping Air Gaps have been corrected")
        else:
            air_gaps[i, :] = np.array([position_tag, airgap_position, airgap_h, c_airgap])
            i += 1

# Characteristic lengths (for mesh sizes)
c_core = core_w/10. *s
c_window = window_w/10 *s
c_conductor = window_w/10 *s

# All points with (x, y, z, mesh_accuracy)
p_outer = np.zeros((4, 4))
p_window = np.zeros((4*n_windows, 4))
p_air_gaps = np.zeros((4*n_air_gaps, 4))

# Different points needed for
if core_type == "EI":

    if y_symmetric == 0:
        # Outer
        p_outer[0][:] = [-(core_w + window_w), -(window_h / 2 + core_w), 0, c_core]
        p_outer[1][:] = [core_w + window_w, -(window_h / 2 + core_w), 0, c_core]
        p_outer[2][:] = [-(core_w + window_w), (window_h / 2 + core_w), 0, c_core]
        p_outer[3][:] = [core_w + window_w, (window_h / 2 + core_w), 0, c_core]
        # Window
        p_window[0] = [-(core_w/2+window_w), -window_h/2, 0, c_window]
        p_window[1] = [-core_w/2, -window_h/2, 0, c_window]
        p_window[2] = [-(core_w/2+window_w)/2, window_h, 0, c_window]
        p_window[3] = [-core_w/2, window_h/2, 0, c_window]
        p_window[4] = [core_w/2, -window_h/2, 0, c_window]
        p_window[5] = [(core_w/2+window_w), -window_h/2, 0, c_window]
        p_window[6] = [core_w/2, window_h/2, 0, c_window]
        p_window[7] = [(core_w/2+window_w), window_h/2, 0, c_window]

    if y_symmetric == 1:
        # Fitting the outer radius to ensure surface area
        r_inner = window_w + core_w/2
        r_outer = np.sqrt((core_w/2)**2 + r_inner**2)  # np.sqrt(window_w**2 + window_w * core_w + core_w**2/2)

        # Outer Core
        # (A_zyl=2pi*r*h => h=0.5r=0.25core_w <=> ensure A_zyl=A_core on the tiniest point)
        p_outer[0][:] = [-r_outer, -(window_h / 2 + core_w/4), 0, c_core]
        p_outer[1][:] = [r_outer, -(window_h / 2 + core_w/4), 0, c_core]
        p_outer[2][:] = [-r_outer, (window_h / 2 + core_w/4), 0, c_core]
        p_outer[3][:] = [r_outer, (window_h / 2 + core_w/4), 0, c_core]

        # Window
        p_window[0] = [-r_inner, -window_h/2, 0, c_window]
        p_window[1] = [-core_w/2, -window_h/2, 0, c_window]
        p_window[2] = [-r_inner, window_h/2, 0, c_window]
        p_window[3] = [-core_w/2, window_h/2, 0, c_window]
        p_window[4] = [core_w/2, -window_h/2, 0, c_window]
        p_window[5] = [r_inner, -window_h/2, 0, c_window]
        p_window[6] = [core_w/2, window_h/2, 0, c_window]
        p_window[7] = [r_inner, window_h/2, 0, c_window]

        # Conductors
        if conductor_type == "full":
            # left window
            # p_conductor[4*i+0][:] = [-r_inner + core_cond_isolation, -window_h/2 + core_cond_isolation, 0, c_conductor]
            # p_conductor[4*i+1][:] = [-core_cond_isolation - core_w/2, -window_h/2 + core_cond_isolation, 0, c_conductor]
            # p_conductor[4*i+2][:] = [-r_inner + core_cond_isolation, window_h/2 - core_cond_isolation, 0, c_conductor]
            # p_conductor[4*i+3][:] = [-core_cond_isolation - core_w/2, window_h/2 - core_cond_isolation, 0, c_conductor]
            # right window
            # full window conductor
            p_conductor[0][:] = [core_cond_isolation + core_w/2, -window_h/2 + core_cond_isolation, 0, c_conductor]
            p_conductor[1][:] = [r_inner - core_cond_isolation, -window_h/2 + core_cond_isolation, 0, c_conductor]
            p_conductor[2][:] = [core_cond_isolation + core_w/2, window_h/2 - core_cond_isolation, 0, c_conductor]
            p_conductor[3][:] = [r_inner - core_cond_isolation, window_h/2 - core_cond_isolation, 0, c_conductor]

        if conductor_type == "stacked":
            # stacking from the ground
            p_conductor = np.empty((4*n_conductors, 4))
            for i in range(0, n_conductors):
                # two conductors above
                p_conductor[4*i+0][:] = [core_cond_isolation + core_w/2, (1-i)*core_cond_isolation + i*(-window_h/2 + core_cond_isolation), 0, c_conductor]
                p_conductor[4*i+1][:] = [r_inner - core_cond_isolation, (1-i)*core_cond_isolation + i*(-window_h/2 + core_cond_isolation), 0, c_conductor]
                p_conductor[4*i+2][:] = [core_cond_isolation + core_w/2, -i*core_cond_isolation + (1-i)*(window_h/2 - core_cond_isolation), 0, c_conductor]
                p_conductor[4*i+3][:] = [r_inner - core_cond_isolation, -i*core_cond_isolation + (1-i)*(window_h/2 - core_cond_isolation), 0, c_conductor]

        if conductor_type == "foil":
            p_conductor = np.empty((4*n_conductors, 4))
            left_bound = core_cond_isolation + core_w/2
            right_bound = r_inner - core_cond_isolation
            x_interpol = np.linspace(left_bound, right_bound, n_conductors+1)
            for i in range(0, n_conductors):
                # Foils
                p_conductor[4 * i + 0][:] = [x_interpol[i] + cond_cond_isolation, -window_h / 2 + core_cond_isolation, 0, c_conductor]
                p_conductor[4 * i + 1][:] = [x_interpol[i+1] - cond_cond_isolation, -window_h / 2 + core_cond_isolation, 0, c_conductor]
                p_conductor[4 * i + 2][:] = [x_interpol[i] + cond_cond_isolation, window_h / 2 - core_cond_isolation, 0, c_conductor]
                p_conductor[4 * i + 3][:] = [x_interpol[i+1] - cond_cond_isolation, window_h / 2 - core_cond_isolation, 0, c_conductor]

    for i in range(0, n_air_gaps):
        # Left leg (-1)
        if air_gaps[i][0] == -1:
            p_air_gaps[i * 4] = [-(core_w + window_w), air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 1] = [-(core_w / 2 + window_w), air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 2] = [-(core_w + window_w), air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 3] = [-(core_w / 2 + window_w), air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]

        # Center leg (0)
        if air_gaps[i][0] == 0:
            p_air_gaps[i * 4] = [-core_w/2, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 1] = [core_w/2, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 2] = [-core_w/2, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 3] = [core_w/2, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]

        # Right leg (+1)
        if air_gaps[i][0] == 1:
            p_air_gaps[i * 4] = [core_w / 2 + window_w, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 1] = [core_w + window_w, air_gaps[i][1] - air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 2] = [core_w / 2 + window_w, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]
            p_air_gaps[i * 4 + 3] = [core_w + window_w, air_gaps[i][1] + air_gaps[i][2] / 2, 0, air_gaps[i][3]]




"""Theoretical inductivity
lz = 0.5
A = lz*(core_w-window_w)/2
mu0 = 1.256*10**(-6)
N = 10
L_theo = None
if n_air_gaps > 0:
    L_theo = N**2 * mu0 * A / airgap_h * 1000
print("Nicht mehr aktuell! Theoretical, idealized inductance: ", L_theo, "mH")
"""