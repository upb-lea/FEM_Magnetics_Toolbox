import numpy as np
import itertools
import matplotlib.pyplot as plt
from femmt_functions import get_dicts_with_keys_and_values, get_dict_with_unique_keys

# Load Nikolas' Data
path = 'C:/Users/tillp/sciebo/Exchange_DAB/current_shapes.npy'
data = np.load(path, allow_pickle=True)
print(data[0])
# N
N1 = [27]  # Turns in main window
N2 = [7]  # Turns in main window
Ns1 = [5]  # Turns in stray window
Ns2 = [6]  # Turns in stray window
N_flat = list(itertools.product(N1, N2, Ns1, Ns2))
N = [np.reshape(N_flat_single, (2, 2)) for N_flat_single in N_flat]

# Goal Parameters
frequency = 250000
L_s1 = 6.8e-5
L_h = 15 * L_s1
n = 3.2
# Inductance Matrix
L_11 = L_s1 + L_h
L_22 = L_h / n**2
M = L_h / n
L = np.array([[L_11, M], [M, L_22]])


# voltage: 175, 235, 295
wp_data = get_dicts_with_keys_and_values(data, power=2500, voltage=235, frequency=frequency,
                                         wp_lambda=L_s1*frequency, wp_n=n, ratio_lm_ls=10)

# print(wp_data)
dict_in = get_dict_with_unique_keys(wp_data, 'wp_ib_il_phase_rad_vec')
dict_m = get_dict_with_unique_keys(wp_data, 'wp_ib_ilm_phase_rad_vec')
dict_out = get_dict_with_unique_keys(wp_data, 'wp_ob_il_phase_rad_vec')


frequencies_in = dict_in['wp_ib_il_frequency_vec']
currents_in = dict_in['wp_ib_il_amplitude_vec']
phases_in = dict_in['wp_ib_il_phase_rad_vec']
time_funct_in = dict_in['wp_ib_il_vec']

time_funct_m = dict_m['wp_ib_il_vec']

frequencies_out = dict_out['wp_ob_il_frequency_vec']
currents_out = dict_out['wp_ob_il_amplitude_vec']
phases_out = dict_out['wp_ob_il_phase_rad_vec']
time_funct_out = dict_out['wp_ob_il_vec']


print(f"{frequencies_in=}\n"
      f"{currents_in=}\n"
      f"{phases_in=}\n"
      f"{time_funct_in=}\n"
      f"{frequencies_out=}\n"
      f"{currents_out=}\n"
      f"{phases_out=}\n"
      f"{time_funct_out=}\n")


# Plotting
figure, axis = plt.subplots(2)

f_switch = frequencies_in[0] if frequencies_in[0] != 0 else frequencies_in[1]
phase_in_1st = phases_in[0] if frequencies_in[0] != 0 else phases_in[1]
phase_out_1st = phases_out[0] if frequencies_out[0] != 0 else phases_out[1]

t = np.linspace(0, 2 * np.pi, 1000) / f_switch
f = frequencies_in

# Time functions
time_funct = np.array(time_funct_in[0])/frequency

current_in = time_funct_in[1]
axis[0].plot(time_funct, current_in, label=f"Time Input")

current_out = np.array(time_funct_out[1])/(-3.2)
axis[0].plot(time_funct, current_out, label=f"Time Output (prim. transformed)")

current_main = current_in + current_out
axis[0].plot(time_funct, current_main, label=f"Time Main")

axis[0].plot(time_funct, time_funct_m[1], label=f"Time Main native")
# ax.plot(t_in, time_funct_in[1], label=f"Main")


# IFFT Output
# Imax_in = currents_in[0]
Imax_in = max(time_funct_in[1])
IFFT_in = 0
for anum, I in enumerate(currents_in):
    phase = phases_in[anum]
    f = frequencies_in[anum]
    IFFT_in += I*np.cos(t*f + phase)

IFFT_in_1st = Imax_in*np.cos(t*f_switch + phase_in_1st)
axis[1].plot(t, IFFT_in_1st, label=f"{f_switch, phase_in_1st}")
# ax.plot(t, I*np.cos(t*f + phase), label=f"{f, phase}")

axis[1].plot(t, IFFT_in, label=f"IFFT Input")

# IFFT Input
# Imax_out = currents_out[0]/3.2
Imax_out = max(time_funct_out[1])/3.2
IFFT_out = 0
for anum, I in enumerate(currents_out):
    phase = phases_out[anum]
    f = frequencies_out[anum]
    IFFT_out += I*np.cos(t*f + phase+np.pi)

IFFT_out_1st = Imax_out*np.cos(t*f_switch + phase_out_1st+np.pi)
axis[1].plot(t, IFFT_out_1st, label=f"{f_switch, phase_out_1st}")
# ax.plot(t, I*np.cos(t*f + phase), label=f"{f, phase}")

axis[1].plot(t, IFFT_out/n, label=f"IFFT Output (prim. transformed)")


# Difference [stray flux]
IFFT_stray = IFFT_in + IFFT_out/n
IFFT_stray_1st = IFFT_in_1st + IFFT_out_1st

axis[1].plot(t, IFFT_stray, label=f"main")
axis[1].plot(t, IFFT_stray_1st, label=f"main_1st")

axis[0].set_ylabel("Current in A")
axis[0].set_xlabel("time in s")
axis[1].set_ylabel("Current in A")
axis[1].set_xlabel("time in s")
axis[0].legend()
axis[1].legend()
plt.show()



I1 = currents_in
I2 = currents_out
I = [I1, I2]

[Phi_top, Phi_bot] = np.matmul(np.matmul(np.linalg.inv(np.transpose(N[0])), L),
                               np.transpose(I))

Phi_stray = np.abs(Phi_top - Phi_bot)

# print(f"{Phi_top=}\n"
#       f"{Phi_bot=}\n"
#       f"{Phi_stray=}\n")

n = -I2/I1


ax = plt.subplot()

plt.plot(n, np.abs(Phi_top), label=f"|Phi_top|")
plt.plot(n, np.abs(Phi_bot), label=f"|Phi_bot|")
plt.plot(n, np.abs(Phi_stray), label=f"|Phi_stray|")

plt.xlabel("n = -I2/I1")
plt.ylabel("Absolute value of magnetic flux in Wb")
plt.legend()
plt.show()
