import numpy as np
from functions import NbrStrands

# -- Control Flags --
flag_imposed_reduced_frequency = 0  # if == 0 --> impose frequency f
flag_excitation_type = 'current'  # 'current', 'current_density', 'voltage'
flag_non_linear_core = 0


# -- Geometry --
core_type = "EI"  # Basic shape of magnetic conductor
#conductor_type = "stacked"  # Vertical packing of conductors
#conductor_type = "full"  # One massive Conductor in each window
#conductor_type = "foil"  # Horizontal packing of conductors
conductor_type = "litz"  # Horizontal packing of conductors
#conductor_type = "solid"  # Horizontal packing of conductors
#conductor_type = None
p_conductor = np.empty(0)  # Case: no conductors only theorretical
y_symmetric = 1  # Mirror-symmetry across y-axis
axi_symmetric = 1  # Axial-symmetric model (idealized full-cylindrical)
n_conductors = 33  # Number of (homogenised) conductors in one window
n_air_gaps = 1  # Number of air gaps - Needs a re-visit for non axi-symmetric case
s = 0.2  # Parameter for mesh-accuracy
# Dimensions in meter
core_cond_isolation = 0.001  # gap between Core and Conductor
cond_cond_isolation = 0.0002  # gap between Core and Conductor
core_w = 0.02  # Axi symmetric case | core_w := core radius
window_w = 0.01  # Winding window width
window_h = 0.03  # Winding window height
if conductor_type == 'solid':
    conductor_radius = 0.0011879

# Litz Approximation
FF = 0.9  # hexagonal packing: ~90.7% are theoretical maximum
n_layers = 6
n_strands = NbrStrands(n_layers)
strand_radius = 0.1e-3
if conductor_type == 'litz':
    conductor_radius = np.sqrt(n_strands/FF)*strand_radius  # Must be calculated from strand
A_cell = np.pi * conductor_radius**2  #* FF  # Surface of the litz approximated hexagonal cell



# -- Materials --
# frequency = 0: mu_rel only used if flag_non_linear_core == 0
# frequency > 0: mu_rel is used
mu0 = 4e-7*np.pi
mu_rel = 3000   # relative Core Permeability
core_material = 95  # 95 := TDK-N95 | Currently only works with Numbers corresponding to BH.pro
sigma = 6e7

# -- Excitation --
# Imposed current, current density or voltage
if flag_excitation_type == 'current':
    current = 3
if flag_excitation_type == 'current_density':
    raise NotImplementedError
if flag_excitation_type == 'voltage':
    voltage = 2
# Frequency and reduced Frequency
frequency = 50000  # in Hz
if flag_imposed_reduced_frequency == 1:
    red_freq = 4
else:
    if frequency != 0:
        delta = np.sqrt(2 / (2 * frequency * np.pi * sigma * mu0))
        red_freq = strand_radius / delta
    else:
        delta = 1e20  # random huge value
        red_freq = 0


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
        airgap_h = 0.0005  #np.random.rand(1)*0.005 + 0.001
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
            x_interpol = np.linspace(left_bound, right_bound, n_conductors+1)  # instead should FF window from inside to outside with fixed copper thickness
            for i in range(0, n_conductors):
                # Foils
                p_conductor[4 * i + 0][:] = [x_interpol[i] + cond_cond_isolation, -window_h / 2 + core_cond_isolation, 0, c_conductor]
                p_conductor[4 * i + 1][:] = [x_interpol[i+1] - cond_cond_isolation, -window_h / 2 + core_cond_isolation, 0, c_conductor]
                p_conductor[4 * i + 2][:] = [x_interpol[i] + cond_cond_isolation, window_h / 2 - core_cond_isolation, 0, c_conductor]
                p_conductor[4 * i + 3][:] = [x_interpol[i+1] - cond_cond_isolation, window_h / 2 - core_cond_isolation, 0, c_conductor]

        if conductor_type == "litz" or conductor_type == "solid":
            p_conductor = []  # center points are stored
            left_bound = core_cond_isolation + core_w/2
            right_bound = r_inner - core_cond_isolation
            top_bound = window_h/2
            bot_bound = -window_h/2

            y = bot_bound+core_cond_isolation+conductor_radius
            x = left_bound + core_cond_isolation + conductor_radius
            i = 0
            # Case n_conductors higher that "allowed" is missing
            while y < top_bound:
                while x < right_bound and i < n_conductors:
                    p_conductor.append([x, y, 0, c_conductor])
                    p_conductor.append([x-conductor_radius, y, 0, c_conductor])
                    p_conductor.append([x, y+conductor_radius, 0, c_conductor])
                    p_conductor.append([x+conductor_radius, y, 0, c_conductor])
                    p_conductor.append([x, y-conductor_radius, 0, c_conductor])
                    i += 1
                    x += conductor_radius * 2 + cond_cond_isolation
                y += conductor_radius * 2 + cond_cond_isolation
                x = left_bound + core_cond_isolation + conductor_radius
            p_conductor = np.asarray(p_conductor)
            if int(p_conductor.shape[0]/5) < n_conductors:
                print("Could not resolve all conductors")

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