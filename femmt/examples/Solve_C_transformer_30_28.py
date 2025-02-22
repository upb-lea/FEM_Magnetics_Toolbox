import numpy as np
import matplotlib.pyplot as plt

# Define the M matrix based on your corrected scenarios

M_case_1 = np.array([
    [1, 0, 0, 1, 0, 1, 1, 0, 0, 0],  # Scenario 1
    [0, 1, 0, -1, 1, 0, 0, 0, 1, 0],  # Scenario 2
    [0, 0, 1, -1, 1, -1, 0, 0, 1, 1],  # Scenario 3
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],  # Scenario 4
    [1, 1, 0, 0, 1, 1, 1, 0, 1, 0],  # Scenario 5
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],  # Scenario 6
    [1, 0, 0, 1, 0, 1, 2, 1, 1, 1],  # Scenario 7
    [0, 1, 1, -2, 2, -1, 0, 0, 2, 1],  # Scenario 8
    [0, 1, 0, -1, 1, 0, 1, 1, 2, 1],  # Scenario 9
    [0, 0, 1, -1, -1, 1, 1, 1, 2, 2],  # Scenario 10
])
M_case_2 = np.array([
    [1, 0, 0, 1, 0, 1, 1, 0, 0, 0],  # Scenario 1
    [-1, 0, -1, 0, -1, 0, 0, 1, 0, 0],  # Scenario 2
    [0, 1, 0, -1, 1, 0, 0, 0, 1, 0],  # Scenario 3
    [0, -1, 1, 0, 0, -1, 0, 0, 0, 1],  # Scenario 4
    [0, 0, -1, 1, -1, 1, -1, 1, 0, 0],  # Scenario 5
    [1, 1, 0, 0, 1, 1, 1, 0, 1, 0],  # Scenario 6
    [1, -1, 1, 1, 0, 0, 1, 0, 0, 1],  # Scenario 7
    [-1, 1, -1, -1, 0, 0, 0, 1, 1, 0],  # Scenario 8
    [-1, -1, 0, 0, -1, -1, 0, 1, 0, 1],  # Scenario 9
    [0, 0, 1, -1, 1, -1, 0, 0, 1, 1],  # Scenario 10
])
# M_case_2 = np.array([
#     [1, 0, 0, 1, 0, 1, 1, 0, 0, 0],  # Scenario 1
#     [1, 0, 1, 0, 1, 0, 0, 1, 0, 0],  # Scenario 2
#     [0, 1, 0, 1, 1, 0, 0, 0, 1, 0],  # Scenario 3
#     [0, 1, 1, 0, 0, 1, 0, 0, 0, 1],  # Scenario 4
#     [0, 0, 1, 1, 1, 1, 1, 1, 0, 0],  # Scenario 5
#     [1, 1, 0, 0, 1, 1, 1, 0, 1, 0],  # Scenario 6
#     [1, 1, 1, 1, 0, 0, 1, 0, 0, 1],  # Scenario 7
#     [1, 1, 1, 1, 0, 0, 0, 1, 1, 0],  # Scenario 8
#     [1, 1, 0, 0, 1, 1, 0, 1, 0, 1],  # Scenario 9
#     [0, 0, 1, 1, 1, 1, 0, 0, 1, 1],  # Scenario 10
#     # [1, 0, 0, 1, 0, 1, 0, 1, 1, 1],
# ])
# M_transpose = M_case_1.T

M_squared1 = M_case_1 ** 2
M_squared2 = M_case_2 ** 2
print("Element-wise Squared Matrix case1:\n", M_squared1)
print("Element-wise Squared Matrix case2:\n", M_squared2)

We_case_1 = np.array([
    2.002418828706215e-11,  # Scenario 1
    1.611847650440696e-11,   # Scenario 2
    4.54639734788541e-11,   # Scenario 3
    1.38508462369603e-11,  # Scenario 4
    5.750010130762671e-12,  # Scenario 5
    1.942906256990097e-11,  # Scenario 6
    4.555462346684122e-11,  # Scenario 7
    1.076387631951943e-10,  # Scenario 8
    3.291112743307101e-11,  # Scenario 9
    6.37630655520892e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])

We_case_2 = np.array([
    2.002418828706215e-11,  # Scenario 1
    1.71520740277683e-11,   # Scenario 2
    1.611847650440696e-11,   # Scenario 3
    1.552613677132804e-11,  # Scenario 4
    5.486657387954085e-11,  # Scenario 5
    5.750010130762671e-12,  # Scenario 6
    1.988388052308087e-11,  # Scenario 7
    2.054869667265151e-11,  # Scenario 8
    4.979462733201061e-12,  # Scenario 9
    4.54639734788541e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])

We_case1 = We_case_1.reshape((10, 1))
We_case2 = We_case_2.reshape((10, 1))
det_M1 = np.linalg.det(M_squared1)
print(f"Determinant of M: {det_M1}")

if det_M1 == 0:
    print("Matrix M is singular and cannot be inverted. Consider revising your simulation scenarios.")
else:
    C_case1 = 2 * np.linalg.inv(M_squared1).dot(We_case1)
    C_case2 = 2 * np.linalg.inv(M_squared2).dot(We_case2)

    # Flatten C to a 1D array for easier indexing
    C_case1 = C_case1.flatten()
    C_case2 = C_case2.flatten()

    # Display
    capacitance_pairs = ["A-B", "C-D", "B-D", "A-C", "B-C",
                         "A-D", "A-E", "B-E", "C-E", "D-E"]
    print("\nIndividual Capacitances:")
    for idx, capacitance in enumerate(C_case1, start=1):
        print(f"C{idx} ({capacitance_pairs[idx-1]}): {capacitance:.5e} F")
    for idx, capacitance in enumerate(C_case2, start=1):
        print(f"C{idx} ({capacitance_pairs[idx-1]}): {capacitance:.5e} F")

    # Calculate the sums for connections in measurment
    C_ABvsCDE = C_case1[2] + C_case1[3] + C_case1[4] + C_case1[5] + C_case1[6] + C_case1[7]
    C_ABEvsCd = C_case1[2] + C_case1[3] + C_case1[4] + C_case1[5] + C_case1[8] + C_case1[9]
    C_ABCDvsE = C_case1[6] + C_case1[7] + C_case1[8] + C_case1[9]
    C_AvsBCDE = C_case1[0] + C_case1[3] + C_case1[5] + C_case1[6]
    C_BvsACDE = C_case1[0] + C_case1[2] + C_case1[4] + C_case1[7]
    C_CvsABDE = C_case1[1] + C_case1[3] + C_case1[4] + C_case1[8]
    C_DvsABCE = C_case1[1] + C_case1[2] + C_case1[5] + C_case1[9]
    C_ACvsBDE = C_case1[0] + C_case1[1] + C_case1[4] + C_case1[5] + C_case1[6] + C_case1[8]
    C_ADvsBCE = C_case1[0] + C_case1[1] + C_case1[2] + C_case1[3] + C_case1[6] + C_case1[9]
    C_BC_ADE = C_case1[0] + C_case1[1] + C_case1[2] + C_case1[3] + C_case1[7] + C_case1[8]

    # Calculate the sums for connections in measurment
    C_ABvsCDE_2 = C_case2[2] + C_case2[3] + C_case2[4] + C_case2[5] + C_case2[6] + C_case2[7]
    C_ABEvsCd_2 = C_case2[2] + C_case2[3] + C_case2[4] + C_case2[5] + C_case2[8] + C_case2[9]
    C_ABCDvsE_2 = C_case2[6] + C_case2[7] + C_case2[8] + C_case2[9]
    C_AvsBCDE_2 = C_case2[0] + C_case2[3] + C_case2[5] + C_case2[6]
    C_BvsACDE_2 = C_case2[0] + C_case2[2] + C_case2[4] + C_case2[7]
    C_CvsABDE_2 = C_case2[1] + C_case2[3] + C_case2[4] + C_case2[8]
    C_DvsABCE_2 = C_case2[1] + C_case2[2] + C_case2[5] + C_case2[9]
    C_ACvsBDE_2 = C_case2[0] + C_case2[1] + C_case2[4] + C_case2[5] + C_case2[6] + C_case2[8]
    C_ADvsBCE_2 = C_case2[0] + C_case2[1] + C_case2[2] + C_case2[3] + C_case2[6] + C_case2[9]
    C_BC_ADE_2 = C_case2[0] + C_case2[1] + C_case2[2] + C_case2[3] + C_case2[7] + C_case2[8]

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

    # Display the sums
    print("\nCalculated Sums:")
    print(f"C_ABvsCDE (C3 + C4 + C5 + C6 + C7 + C8): {C_ABvsCDE_2:.5e} F")
    print(f"C_ABEvsCd (C3 + C4 + C5 + C6 + C9 + C10): {C_ABEvsCd_2:.5e} F")
    print(f"C_ABCDvsE (C7 + C8 + C9 + C10): {C_ABCDvsE_2:.5e} F")
    print(f"C_AvsBCDE (C1 + C4 + C6 + C7): {C_AvsBCDE_2:.5e} F")
    print(f"C_BvsACDE (C1 + C3 + C5 + C8): {C_BvsACDE_2:.5e} F")
    print(f"C_CvsABDE (C2 + C4 + C5 + C9): {C_CvsABDE_2:.5e} F")
    print(f"C_DvsABCE (C2 + C3 + C6 + C10): {C_DvsABCE_2:.5e} F")
    print(f"C_ACvsBDE (C1 + C2 + C5 + C6 + C7 + C9): {C_ACvsBDE_2:.5e} F")
    print(f"C_ADvsBCE (C1 + C2 + C3 + C4 + C7 + C10): {C_ADvsBCE_2:.5e} F")
    print(f"C_BC_ADE (C1 + C2 + C3 + C4 + C8 + C9): {C_BC_ADE_2:.5e} F")

# new
energy_data = [We_case1]
energy_labels = ["Case1"]

# Calculate C for each case
capacitances = []
for We in energy_data:
    We = We.reshape((10, 1))
    M_squared = M_case_1 ** 2
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
# connection_sums = {
#     'C_ABvsCDE': lambda C: C[2] + C[3] + C[4] + C[5] + C[6] + C[7],
#     'C_ABEvsCd': lambda C: C[2] + C[3] + C[4] + C[5] + C[8] + C[9],
#     'C_ABCDvsE': lambda C: C[6] + C[7] + C[8] + C[9],
#     'C_AvsBCDE': lambda C: C[0] + C[3] + C[5] + C[6],
#     'C_BvsACDE': lambda C: C[0] + C[2] + C[4] + C[7],
#     'C_CvsABDE': lambda C: C[1] + C[3] + C[4] + C[8],
#     'C_DvsABCE': lambda C: C[1] + C[2] + C[5] + C[9],
#     'C_ACvsBDE': lambda C: C[0] + C[1] + C[4] + C[5] + C[6] + C[8],
#     'C_ADvsBCE': lambda C: C[0] + C[1] + C[2] + C[3] + C[6] + C[9],
#     'C_BC_ADE': lambda C: C[0] + C[1] + C[2] + C[3] + C[7] + C[8],
#     'C_BD_ACE': lambda C: C[0] + C[1] + C[4] + C[5] + C[7] + C[9],
# }
connection_sums = {
    'C_ABvsCDE': lambda C: C[2] + C[3] + C[4] + C[5] + C[6] + C[7],
    'C_ABEvsCd': lambda C: C[2] + C[3] + C[4] + C[5] + C[8] + C[9],
    'C_ABCDvsE': lambda C: C[6] + C[7] + C[8] + C[9],
    'C_AvsBCDE': lambda C: C[0] + C[3] + C[5] + C[6],
    'C_BvsACDE': lambda C: C[0] + C[2] + C[4] + C[7],
    'C_DvsABCE': lambda C: C[1] + C[3] + C[4] + C[8],
    'C_CvsABDE': lambda C: C[1] + C[2] + C[5] + C[9],
    'C_ADvsBCE': lambda C: C[0] + C[1] + C[4] + C[5] + C[6] + C[8],
    'C_ACvsBDE': lambda C: C[0] + C[1] + C[2] + C[3] + C[6] + C[9],
    'C_BD_ACE': lambda C: C[0] + C[1] + C[2] + C[3] + C[7] + C[8],
    'C_BC_ADE': lambda C: C[0] + C[1] + C[4] + C[5] + C[7] + C[9],
}


connection_measurement = {
    'C_ABvsCDE': 1.51e-10,
    'C_ABEvsCd': 1.35e-10,
    'C_ABCDvsE': 3.21e-11,
    'C_AvsBCDE': 1.8145e-10,
    'C_BvsACDE': 1.9565e-10,
    'C_CvsABDE': 2.0014e-10,
    'C_DvsABCE': 1.4616e-10,
    'C_ACvsBDE': 1.3501e-11,
    'C_ADvsBCE': 1.3021e-10,
    'C_BC_ADE': 9.0531e-10,
    'C_BD_ACE':1.3806e-11,
}

# Plot capacitance values as discrete points
plt.figure(figsize=(12, 6))
markers = ['o', 's', '^', 'D', '*', 'v', 'x']
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
    C = np.linalg.inv(M_case_1 ** 2).dot(2 * We).flatten()
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

condition_number = np.linalg.cond(M_squared)
print(f"Condition number of M_squared: {condition_number}")

# Define a function to calculate, print, and plot connection sums with error
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


# Iterate over all simulations and calculate sums and errors
plt.figure(figsize=(12, 6))
for energy, label in zip(energy_data, energy_labels):
    We = energy.reshape((10, 1))
    C = np.linalg.inv(M_case_1 ** 2).dot(2 * We).flatten()

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

# Iterate over all simulations and calculate sums and errors
plt.figure(figsize=(12, 6))

for energy, label in zip(energy_data, energy_labels):
    We = energy.reshape((10, 1))
    C = np.linalg.inv(M_case_1 ** 2).dot(2 * We).flatten()

    # Print the connection sums, measured values, and errors
    errors = []
    calculated_sums = [connection_sums[key](C) for key in connection_sums]
    measured_sums = list(connection_measurement.values())

    for i, (calc_sum, meas_sum) in enumerate(zip(calculated_sums, measured_sums)):
        error = abs(calc_sum - meas_sum) / calc_sum * 100  # Calculate error percentage
        errors.append(error)
        plt.text(i, calc_sum, f"{error:.1f}%", fontsize=10, ha='center', va='bottom', color='blue')  # Annotate error

    plt.scatter(range(len(calculated_sums)), calculated_sums, label=f"Calculated - {label}", marker='o')

# Plot the measured connections
plt.scatter(range(len(measured_sums)), measured_sums, label="Measured Connections", color="red", marker='x')

# Set plot details
plt.xticks(range(len(connection_sums)), list(connection_sums.keys()), rotation=45)
plt.ylabel('Sum of Capacitances (F)')
plt.title('Connection Sums and Errors for Different Cases')
plt.legend()
plt.grid(True)
plt.show()

