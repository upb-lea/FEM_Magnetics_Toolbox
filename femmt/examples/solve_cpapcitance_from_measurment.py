import numpy as np
import matplotlib.pyplot as plt

# Define the M matrix based on your corrected scenarios
M = np.array([
    [1, 0, 0, 1, 0, 1, 1, 0, 0, 0],  # Scenario 1
    [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],  # Scenario 2
    [0, 1, 0, 1, 1, 0, 0, 0, 1, 0],  # Scenario 3
    [0, 1, 1, 0, 0, 1, 0, 0, 0, 1],  # Scenario 4
    [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],  # Scenario 5
    [1, 1, 0, 0, 1, 1, 1, 0, 1, 0],  # Scenario 6
    [1, 1, 1, 1, 0, 0, 1, 0, 0, 1],  # Scenario 7
    [1, 1, 1, 1, 0, 0, 0, 1, 1, 0],  # Scenario 8
    [1, 1, 0, 0, 1, 1, 0, 1, 0, 1],  # Scenario 9
    [0, 0, 1, 1, 1, 1, 0, 0, 1, 1],  # Scenario 10
    # [1, 0, 0, 1, 0, 1, 0, 1, 1, 1],
])

M_squared = M ** 2
print("Element-wise Squared Matrix M:\n", M_squared)
# Define the energy vector
# with error 50%
# adjusting
We_no_bobbin = np.array([
    1.107533022806718e-11,  # Scenario 1
    2.998536727832599e-11,   # Scenario 2
    2.441518883756004e-11,   # Scenario 3
    1.126487136125143e-11,  # Scenario 4
    3.436564097908541e-11,  # Scenario 5
    3.11986239282366e-11,  # Scenario 6
    2.014172856588652e-11,  # Scenario 7
    1.534054057620717e-11,  # Scenario 8
    2.924182091720404e-11,  # Scenario 9
    3.222388455288184e-11,   # Scenario 10
    # 1.114304711308078e-11
])
#with bobbin in femmt
We_with_bobbin = np.array([
    2.358792372177351e-11,  # Scenario 1
    3.10183406267461e-11,   # Scenario 2
    2.56393791428112e-11,   # Scenario 3
    1.248543827275428e-11,  # Scenario 4
    5.172790534894289e-11,  # Scenario 5
    4.484317223246496e-11,  # Scenario 6
    3.387134997216844e-11,  # Scenario 7
    1.752841658974935e-11,  # Scenario 8
    3.15035211362564e-11,  # Scenario 9
    3.478179986637866e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])
# with modified boobin
We_with_modified_bobbin = np.array([
    2.28445470419967E-11,  # Scenario 1
    3.052061386018125e-11,   # Scenario 2
    2.462827943671277e-11,   # Scenario 3
    1.178675793204942e-11,  # Scenario 4
    5.01055464536329e-11,  # Scenario 5
    4.315723991354992e-11,  # Scenario 6
    3.241385446834553e-11,  # Scenario 7
    1.614028079446477e-11,  # Scenario 8
    3.028421087334619e-11,  # Scenario 9
    3.323122435810969e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])
# with modified left
We_with_modified_bobbin_left = np.array([
    1.245261167911289e-11,  # Scenario 1
    3.342773662349677e-11,   # Scenario 2
    3.147191424565136e-11,   # Scenario 3
    1.775552140732446e-11,  # Scenario 4
    3.858201796541978e-11,  # Scenario 5
    3.8979854718569e-11,  # Scenario 6
    2.770154315880028e-11,  # Scenario 7
    2.210290692434063e-11,  # Scenario 8
    3.793604217470713e-11,  # Scenario 9
    3.827079696286146e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])
# with_modified_permitivity_air
We_with_modified_permitivity_air = np.array([
    2.82304511778657E-11,  # Scenario 1
    8.25872891113608e-11,   # Scenario 2
    7.706769180919169e-11,   # Scenario 3
    4.266150024100241e-11,  # Scenario 4
    9.169464496063301e-11,  # Scenario 5
    9.339105676653725e-11,  # Scenario 6
    6.480989614325092e-11,  # Scenario 7
    5.295661549437068e-11,  # Scenario 8
    9.231609431861514e-11,  # Scenario 9
    9.224172711066144e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])
# with_modified_new_scheme_cond
We_with_modified_new_scheme_cond = np.array([
    1.509699848504535e-11,  # Scenario 1
    4.81099275325166E-11,   # Scenario 2
    4.264892437955988e-11,   # Scenario 3
    2.043270191945802e-11,  # Scenario 4
    5.421253505271627e-11,  # Scenario 5
    5.022384218284499e-11,  # Scenario 6
    3.183455552006216e-11,  # Scenario 7
    2.718229984170319e-11,  # Scenario 8
    4.790872894200632e-11,  # Scenario 9
    5.27967519750478e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])
We_with_modified_new_scheme_cond_with_ins_delta = np.array([
    1.601968528078265e-11,  # Scenario 1
    3.965640394258638e-11,   # Scenario 2
    4.090717748855067e-11,   # Scenario 3
    2.393673718687787e-11,  # Scenario 4
    4.247113724840581e-11,  # Scenario 5
    5.122712780820113e-11,  # Scenario 6
    3.72302024153097e-11,  # Scenario 7
    3.252496641928966e-11,  # Scenario 8
    4.912307419319937e-11,  # Scenario 9
    4.086817229316922e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])
We_Othman = np.array([
    2.003333780181938e-11,  # Scenario 1
    3.1950089708645e-11,   # Scenario 2
    3.412096497301497e-11,   # Scenario 3
    2.106331013737882e-11,  # Scenario 4
    4.120287404294935e-11,  # Scenario 5
    4.976356311074315e-11,  # Scenario 6
    3.899419055671076e-11,  # Scenario 7
    2.932421372288165e-11,  # Scenario 8
    4.189076325883805e-11,  # Scenario 9
    3.15589983085852e-11,   # Scenario 10
    # 2.350977482015328e-11    # 11
])

# Ensure We is a column vector
# We = We_no_bobbin.reshape((10, 1))
# We = We_with_bobbin.reshape((10, 1))
# We = We_with_modified_bobbin.reshape((10, 1))
We = We_Othman.reshape((10, 1))
# Check if M is invertible by computing its determinant
det_M = np.linalg.det(M_squared)
print(f"Determinant of M: {det_M}")

if det_M == 0:
    print("Matrix M is singular and cannot be inverted. Consider revising your simulation scenarios.")
else:
    # Solve for C: C = M^{-1} * (2 * We)
    C = np.linalg.inv(M_squared).dot(2 * We)

    # Flatten C to a 1D array for easier indexing
    C = C.flatten()

    # Display
    capacitance_pairs = ["A-B", "C-D", "A-C", "B-D", "B-C",
                         "A-D", "A-E", "B-E", "C-E", "D-E"]
    print("\nIndividual Capacitances:")
    for idx, capacitance in enumerate(C, start=1):
        print(f"C{idx} ({capacitance_pairs[idx-1]}): {capacitance:.5e} F")

    # Calculate the sums for connections in measurment
    C_ABvsCDE = C[2] + C[3] + C[4] + C[5] + C[6] + C[7]
    C_ABEvsCd = C[2] + C[3] + C[4] + C[5] + C[8] + C[9]
    C_ABCDvsE = C[6] + C[7] + C[8] + C[9]
    C_AvsBCDE = C[0] + C[3] + C[5] + C[6]
    C_BvsACDE = C[0] + C[2] + C[4] + C[7]
    C_CvsABDE = C[1] + C[3] + C[4] + C[8]
    C_DvsABCE = C[1] + C[2] + C[5] + C[9]
    C_ACvsBDE = C[0] + C[1] + C[4] + C[5] + C[6] + C[8]
    C_ADvsBCE = C[0] + C[1] + C[2] + C[3] + C[6] + C[9]
    C_BC_ADE = C[0] + C[1] + C[2] + C[3] + C[7] + C[8]

    # Display the sums
    print("\nCalculated Sums:")
    print(f"C_ABvsCDE (C3 + C4 + C5 + C6 + C7 + C8): {C_ABvsCDE:.5e} F")
    print(f"C_ABEvsCd (C3 + C4 + C5 + C6 + C9 + C10): {C_ABEvsCd:.5e} F")
    print(f"C_ABCDvsE (C7 + C8 + C9 + C10): {C_ABCDvsE:.5e} F")
    print(f"C_AvsBCDE (C1 + C4 + C6 + C7): {C_AvsBCDE:.5e} F")
    print(f"C_BvsACDE (C1 + C3 + C5 + C8): {C_BvsACDE:.5e} F")
    print(f"C_CvsABDE (C2 + C4 + C5 + C9): {C_CvsABDE:.5e} F")
    print(f"C_DvsABCE (C2 + C3 + C6 + C10): {C_DvsABCE:.5e} F")
    print(f"C_ACvsBDE (C1 + C2 + C5 + C6 + C7 + C9): {C_ACvsBDE:.5e} F")
    print(f"C_ADvsBCE (C1 + C2 + C3 + C4 + C7 + C10): {C_ADvsBCE:.5e} F")
    print(f"C_BC_ADE (C1 + C2 + C3 + C4 + C8 + C9): {C_BC_ADE:.5e} F")

# new
energy_data = [We_no_bobbin, We_with_bobbin, We_with_modified_bobbin, We_with_modified_bobbin_left, We_with_modified_permitivity_air, We_with_modified_new_scheme_cond, We_with_modified_new_scheme_cond_with_ins_delta, We_Othman]
energy_labels = ["No Bobbin", "With Bobbin", "Modified Bobbin Top and Bot", "Modified Bobbin Top, Bot and Left", "Modified_per_of_air", "modified_new_scheme_cond", "modified_new_scheme_cond_with_ins_delta", "We_Othman"]

# Calculate C for each case
capacitances = []
for We in energy_data:
    We = We.reshape((10, 1))
    M_squared = M ** 2
    C = np.linalg.inv(M_squared).dot(2 * We).flatten()
    capacitances.append(C)
print(capacitances)
# Plot the energy values
plt.figure(figsize=(12, 6))
for idx, We in enumerate(energy_data):
    plt.plot(range(1, 11), We, label=energy_labels[idx])
plt.xlabel("Scenario")
plt.ylabel("Energy (J)")
plt.title("Energy Values for Different Cases")
plt.legend()
plt.grid(True)
plt.show()



# Plot the calculated sums
connection_sums = {
    'C_ABvsCDE': lambda C: C[2] + C[3] + C[4] + C[5] + C[6] + C[7],
    'C_ABEvsCd': lambda C: C[2] + C[3] + C[4] + C[5] + C[8] + C[9],
    'C_ABCDvsE': lambda C: C[6] + C[7] + C[8] + C[9],
    'C_AvsBCDE': lambda C: C[0] + C[3] + C[5] + C[6],
    'C_BvsACDE': lambda C: C[0] + C[2] + C[4] + C[7],
    'C_CvsABDE': lambda C: C[1] + C[3] + C[4] + C[8],
    'C_DvsABCE': lambda C: C[1] + C[2] + C[5] + C[9],
    'C_ACvsBDE': lambda C: C[0] + C[1] + C[4] + C[5] + C[6] + C[8],
    'C_ADvsBCE': lambda C: C[0] + C[1] + C[2] + C[3] + C[6] + C[9],
    'C_BC_ADE': lambda C: C[0] + C[1] + C[2] + C[3] + C[7] + C[8],
    'C_BD_ACE': lambda C: C[0] + C[1] + C[4] + C[5] + C[7] + C[9],
}
connection_measurement = {
    'C_ABvsCDE': 1.41e-10,
    'C_ABEvsCd': 1.23e-10,
    'C_ABCDvsE': 3.81E-11,
    'C_AvsBCDE': 2.86E-11,
    'C_BvsACDE': 7.20E-11,
    'C_CvsABDE': 3.95E-11,
    'C_DvsABCE': 3.87E-11,
    'C_ACvsBDE': 8.13E-11,
    'C_ADvsBCE': 6.58E-11,
    'C_BC_ADE': 6.59E-11,
    'C_BD_ACE':8.11E-11,
}


# Plot capacitance values as discrete points
plt.figure(figsize=(12, 6))
markers = ['o', 's', '^', 'D', '*', 'v', 'x','p']
for idx, C in enumerate(capacitances):
    plt.scatter(range(1, 11), C, label=f"{energy_labels[idx]}", marker=markers[idx])
    # Annotate each capacitance value
    for i, c_value in enumerate(C):
        plt.text(i + 1, c_value, f"C{i+1}", fontsize=8, ha='right')

plt.xlabel("Capacitance Index")
plt.ylabel("Capacitance (F)")
plt.title("Capacitance Values for Different Cases (Discrete Points)")
plt.legend()
plt.grid(True)
plt.show()


# Plot the calculated sums for each energy case
for energy, label in zip(energy_data, energy_labels):
    We = energy.reshape((10, 1))
    C = np.linalg.inv(M ** 2).dot(2 * We).flatten()
    sums = [connection_sums[key](C) for key in connection_sums]
    plt.scatter(list(range(len(connection_sums))), sums, label=f"Calculated - {label}")
    # for i, sum_value in enumerate(sums):
    #     plt.text(i, sum_value, f"{list(connection_sums.keys())[i]}", fontsize=8, ha='center')

# Plot the measured connections
measured_sums = list(connection_measurement.values())
plt.scatter(list(range(len(connection_measurement))), measured_sums, label="Measured Connections", color="red", marker='x')
# for i, measured_value in enumerate(measured_sums):
#     plt.text(i, measured_value, f"{list(connection_measurement.keys())[i]}", fontsize=8, color="red", ha='center')

# Set plot details
plt.xticks(range(len(connection_sums)), list(connection_sums.keys()), rotation=45)
plt.ylabel('Sum of Capacitances (F)')
plt.title('Connection Sums for Different Cases and Measured Connections')
plt.legend()
plt.grid(True)
plt.show()

def print_and_plot_sums_with_error(C, label):
    print(f"\n--- Connection Sums for Simulation: {label} ---")
    errors = []

    for name, func in connection_sums.items():
        sum_value = func(C)
        measured_value = connection_measurement[name]
        error = abs(sum_value - measured_value) / sum_value * 100
        errors.append(error)

        print(f"{name}: {sum_value:.5e} F (Measured: {measured_value:.5e} F) | Error: {error:.2f}%")

    return errors

#Iterate over all simulations and calculate sums and errors
plt.figure(figsize=(12, 6))
for energy, label in zip(energy_data, energy_labels):
    We = energy.reshape((10, 1))
    C = np.linalg.inv(M ** 2).dot(2 * We).flatten()

    # Print the connection sums, measured values, and errors
    errors = print_and_plot_sums_with_error(C, label)

    # Plot the calculated and measured sums
    calculated_sums = [connection_sums[key](C) for key in connection_sums]
    measured_sums = list(connection_measurement.values())

    plt.scatter(range(len(calculated_sums)), calculated_sums, label=f"Calculated - {label}", marker='o')
    plt.scatter(range(len(measured_sums)), measured_sums, label="Measured Connections", color="red", marker='x')

# Set plot details
plt.xticks(range(len(connection_sums)), list(connection_sums.keys()), rotation=45)
plt.ylabel('Sum of Capacitances (F)')
plt.title('Connection Sums and Errors for Different Cases')
plt.legend()
plt.grid(True)
plt.show()

