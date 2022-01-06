from femmt import MagneticComponent
import re
import itertools
import numpy as np

# -- Settings/Flags --
accuracy = 0.5
max_iter = 1

# -- Fixed values --
frequency = 100000
I0 = 10

# -- Bounds --
L_11_min = 9.41e-06
L_11_max = 9.4e-06
L_22_min = 1.9e-05
L_22_max = 2.6e-05
M_min = -1.34e-05
M_max = -1.28e-05

# -- Parameter Vector --
# TODO: randomized Parameter Vectors
N1 = [5, 6, 7]
N2 = [9, 10, 11]
window_h = [0.025, 0.030, 0.035]
window_w = [0.08, 0.01, 0.012]
core_w = [0.018, 0.02, 0.022]
air_gap_h = [0.0019, 0.0020, 0.0021]
conductor_radius1 = [0.0014, 0.0015, 0.0016, 0.0017, 0.0018]
conductor_radius2 = [0.0014, 0.0015, 0.0016, 0.0017, 0.0018]
# x_ = [N1, N2, window_h, window_w, core_w, air_gap_h, conductor_radius1, conductor_radius2]
x__ = list(itertools.product(N1, N2, window_h, window_w, core_w, air_gap_h, conductor_radius1, conductor_radius2))

# -- Results/Logs --
parameters = []
results = []

# -- Brute Force --
count = 0
for x_ in x__:

    # -- Component Selection --
    geo = MagneticComponent(component_type="transformer")

    # -- Core --
    geo.update_core(core_type="EI", window_h=x_[2], window_w=x_[3], core_w=x_[4])

    # -- Air Gaps --
    geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[x_[5]])

    # -- Simulation --
    geo.update_conductors(n_turns=[x_[0], x_[1]],
                          conductor_type=["solid", "solid"],
                          conductor_radix=[x_[6], x_[7]],
                          scheme=["square", "square"],
                          cond_cond_isolation=[0.0002, 0.0002, 0.0008],
                          core_cond_isolation=[0.0004, 0.0004],
                          )

    # Start Simulation
    geo.get_inductances(I0=I0, op_frequency=frequency, skin_mesh_factor=accuracy, visualize=0)

    if geo.valid:
        # if L_11_min < geo.L_11 < L_11_max and L_22_min < geo.L_22 < L_22_max and M_min < geo.M < M_max:
        if True:
            # geo.single_simulation(freq=frequency, current=[I0, N2/N1*I0], phi=[0, 0], skin_mesh_factor=accuracy)
            geo.excitation(f=frequency, i=[I0, x_[1]/x_[0]*I0], phases=[0, 0])  # frequency and current
            geo.file_communication()
            geo.pre_simulate()
            geo.simulate()
            # geo.visualize()
            # learn sth
            data = geo.get_loss_data(loss_type="solid_loss", last_n_values=1)
            for lines in data:
                print(re.findall(r"[-+]?\d*\.\d+|\d+", lines[0]))
                fls = re.findall(r"[-+]?\d*\.\d+|\d+", lines[0])
                Pv = float(fls[1])
            print(f"Power Loss: {Pv}")
            results.append([Pv, geo.L_11, geo.L_22, geo.M])

            parameters.append(x_)

            count += 1
            print(f"Run no. {count}")

            if count % 10 == 0:
                results_tmp = np.asarray(results)
                np.save("results_tmp", results_tmp)
                parameters_tmp = np.asarray(parameters)
                np.save("parameters_tmp", parameters_tmp)
        else:
            print(f"Inductances are not all in bounds! \n {geo.L_11, geo.L_22, geo.M}")


results = np.asarray(results)
np.save("results", results)
parameters = np.asarray(parameters)
np.save("parameters", parameters)

# Zwischendurch abspeichern + Mitzählen

"""
# -- Parameter Vector --
N1 = [6]
N2 = [10]
window_h = [0.03]
window_w = [0.01]
core_w = [0.02]
air_gap_h = [0.002]
conductor_radius1 = [0.0016]
conductor_radius2 = [0.0015, 0.0014]
"""