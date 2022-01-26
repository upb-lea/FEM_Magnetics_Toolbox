import itertools
import numpy as np

# -------------------------------------
# Reluctance Model Parameters
window_h = [0.0295, 12, 575]
window_w = [0.012, 4, 6]
core_w = [0.015]
midpoint_percent = [25, 30, 35]
N1 = [27, 28]  # Turns in main window
N2 = [7]  # Turns in main window
Ns1 = [5]  # Turns in stray window
Ns2 = [6]  # Turns in stray window
N_flat = list(itertools.product(N1, N2, Ns1, Ns2))
print(N_flat)

N = [np.array_repr(np.reshape(N_flat_single, (2, 2))).replace('\n', '') for N_flat_single in N_flat]
# N = [np.reshape(N_flat_single, (2, 2)) for N_flat_single in N_flat]
#N = [[[N1, N2],
#     [Ns1, Ns2]]]
print(N)

# -------------------------------------
# Create List of Dictionaries for Reluctance Model
reluctance_parameter_categories = ["window_h", "window_w", "core_w", "midpoint_percent", "N"]
reluctance_parameter_values = list(itertools.product(window_h, window_w, core_w, midpoint_percent, N))
print(reluctance_parameter_values)

reluctance_parameters = []
for objects in reluctance_parameter_values:
    reluctance_parameters.append({key: value for key, value in zip(reluctance_parameter_categories, objects)})

print(reluctance_parameters)


# -------------------------------------
# This is what happens with the invalid parameter configurations in the Reluctance Model:
# Invalid Parameter Sets are replaced with a None
# Valid Parameter Sets are extended with the air gap lengths [only shown for one arbitrary example]
reluctance_parameters[0] = None
reluctance_parameters[4] = None
reluctance_parameters[1]["some_airgap_length"] = 0.001
#print(reluctance_parameters)
print(len(reluctance_parameters))

# -------------------------------------
# Filter all entries, that are None
# Elements of list reluctance_parameters are either a dictionary or None
valid_reluctance_parameters = [x for x in reluctance_parameters if x is not None]
#print(valid_reluctance_parameters)
print(len(valid_reluctance_parameters))


# -------------------------------------
# Strand Parameters
strand_radius = [0.025e-3]
N_strands_prim = [300, 450]
N_strands_sec = [300, 450]

# -------------------------------------
# Create List of Dictionaries for FEM simulations
non_reluctance_categories = ["strand_radius", "N_strands_prim", "N_strands_sec"]
non_reluctance_values = []
non_reluctance_values = list(itertools.product(strand_radius, N_strands_prim, N_strands_sec))
#print(non_reluctance_values)


non_reluctance_parameters = []
for objects in non_reluctance_values:
    non_reluctance_parameters.append({key: value for key, value in zip(non_reluctance_categories, objects)})

#print(non_reluctance_parameters)
print(len(non_reluctance_parameters))


# -------------------------------------
# Bring together the
FEM_parameters = []
for reluctance_parameters in valid_reluctance_parameters:
    for non_reluctance_parameter in non_reluctance_parameters:
        print(non_reluctance_parameter)
        print(reluctance_parameters)
        FEM_parameters.append(dict(reluctance_parameters, **non_reluctance_parameter))

#print(FEM_parameters)
print(len(FEM_parameters))


# -------------------------------------
# Manually Access/Manipulate the dictionaries
# reluctance_parameters[0]["window_w"] = "something"
# print(reluctance_parameters)




