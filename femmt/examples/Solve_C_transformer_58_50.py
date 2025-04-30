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
    2.430180364410492e-11,  # Scenario 1
    3.963870021343636e-11,   # Scenario 2
    4.628271783338096e-11,   # Scenario 3
    1.705897537969348e-11,  # Scenario 4
    5.084215412258796e-11,  # Scenario 5
    5.061364115385603e-11,  # Scenario 6
    6.041130951808294e-11,  # Scenario 7
    1.521894917175105e-10,  # Scenario 8
    5.994497323888096e-11,  # Scenario 9
    7.127667730542149e-11,   # Scenario 10
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
    # C_BD_ACE = C_case1[0] + C_case1[1] + C_case1[4] + C_case1[5] + C_case1[7] + C_case1[9]

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
    #print(f"C_BD_ACE (C1 + C2 + C5 + C6 + C8 + C10): {C_BC_ADE:.5e} F")

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

#Plot the calculated sums
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
    #'C_BD_ACE': lambda C: C[0] + C[1] + C[4] + C[5] + C[7] + C[9],
}
# connection_sums = {
#     'C_ABvsCDE': lambda C: C[2] + C[3] + C[4] + C[5] + C[6] + C[7],
#     'C_ABEvsCd': lambda C: C[2] + C[3] + C[4] + C[5] + C[8] + C[9],
#     'C_ABCDvsE': lambda C: C[6] + C[7] + C[8] + C[9],
#     'C_AvsBCDE': lambda C: C[0] + C[3] + C[5] + C[6],
#     'C_BvsACDE': lambda C: C[0] + C[2] + C[4] + C[7],
#     'C_DvsABCE': lambda C: C[1] + C[3] + C[4] + C[8],
#     'C_CvsABDE': lambda C: C[1] + C[2] + C[5] + C[9],
#     'C_ADvsBCE': lambda C: C[0] + C[1] + C[4] + C[5] + C[6] + C[8],
#     'C_ACvsBDE': lambda C: C[0] + C[1] + C[2] + C[3] + C[6] + C[9],
#     'C_BD_ACE': lambda C: C[0] + C[1] + C[2] + C[3] + C[7] + C[8],
#     'C_BC_ADE': lambda C: C[0] + C[1] + C[4] + C[5] + C[7] + C[9],
# }


connection_measurement = {
    'C_ABvsCDE': 1.257189992699E-10,
    'C_ABEvsCd': 1.036613057794E-10,
    'C_ABCDvsE': 3.581451055984E-11,
    'C_AvsBCDE': 3.515076277291E-11,
    'C_BvsACDE': 7.638724988931E-11,
    'C_CvsABDE': 7.249036736151E-11,
    'C_DvsABCE': 2.855875033066E-11,
    'C_ACvsBDE': 1.220712831988E-10,
    'C_ADvsBCE': 8.111945395629E-11,
    'C_BC_ADE': 6.237211347932E-11,
    # 'C_BD_ACE':1.3806e-11,
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

condition_number = np.linalg.cond(M_squared1)
print(f"Condition number of M_squared: {condition_number}")

# Define a function to calculate, print, and plot connection sums with error
def print_and_plot_sums_with_error(C, label):
    print(f"\n--- Connection Sums for Simulation: {label} ---")
    errors = []

    for name, func in connection_sums.items():
        sum_value = func(C)
        measured_value = connection_measurement[name]
        error = abs(measured_value - sum_value) / measured_value * 100
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
        error = abs(meas_sum - calc_sum) / meas_sum * 100  # Calculate error percentage
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

# Iterate over all simulations and calculate sums and error ratios
plt.figure(figsize=(12, 6))

for energy, label in zip(energy_data, energy_labels):
    We = energy.reshape((10, 1))
    C = np.linalg.inv(M_case_1 ** 2).dot(2 * We).flatten()

    # Calculate the connection sums for simulation
    calculated_sums = [connection_sums[key](C) for key in connection_sums]
    measured_sums = list(connection_measurement.values())

    # Compute error ratio and annotate each point
    for i, (calc_sum, meas_sum) in enumerate(zip(calculated_sums, measured_sums)):
        # error_ratio = meas_sum / calc_sum  # Calculate error ratio
        error_ratio = calc_sum / meas_sum  # Calculate error ratio
        plt.text(i, calc_sum, f"{error_ratio:.2f}", fontsize=10, ha='center', va='bottom', color='blue')

    # Scatter plot for calculated sums
    plt.scatter(range(len(calculated_sums)), calculated_sums, label=f"Simulated - {label}", marker='o')

# Scatter plot for measured sums
plt.scatter(range(len(measured_sums)), measured_sums, label="Measured Connections", color="red", marker='x')

# Set plot details
plt.xticks(range(len(connection_sums)), list(connection_sums.keys()), rotation=45)
plt.ylabel('Sum of Capacitances (F)')
plt.title('Connection Sums and Error Ratios for Different Cases')
plt.legend()
plt.grid(True)
plt.show()


# Define clean terminal labels for the x-axis
connection_labels_clean = [
    "AB vs CDE", "ABE vs CD", "ABCD vs E", "A vs BCDE", "B vs ACDE",
    "C vs ABDE", "D vs ABCE", "AC vs BDE", "AD vs BCE", "BC vs ADE"
]

# Focused on Case 1: Plot connection sums and annotate error ratios
plt.figure(figsize=(12, 6))

for energy, label in zip(energy_data, energy_labels):
    We = energy.reshape((10, 1))
    C = np.linalg.inv(M_case_1 ** 2).dot(2 * We).flatten()

    # Calculate simulated and measured sums
    calculated_sums = [connection_sums[key](C) for key in connection_sums]
    measured_sums = list(connection_measurement.values())

    # Annotate error ratio (simulated / measured)
    for i, (calc_sum, meas_sum) in enumerate(zip(calculated_sums, measured_sums)):
        error_ratio = calc_sum / meas_sum
        plt.text(i, calc_sum, f"{error_ratio:.2f}×", fontsize=9, ha='center', va='bottom', color='blue')

    # Plot simulated values
    plt.scatter(range(len(calculated_sums)), calculated_sums,
                label="Simulated Capacitance Sum (Case 1)", color='blue', marker='o')

# Plot measured values
plt.scatter(range(len(measured_sums)), measured_sums,
            label="Measured Capacitance Sum", color='red', marker='x')

# Plot formatting
plt.xticks(range(len(connection_labels_clean)), connection_labels_clean, rotation=45)
plt.ylabel('Capacitance Sum (F)')
plt.title('Comparison of Simulated and Measured Capacitance Sums with Error Ratios')
plt.legend(loc='upper right')
plt.tight_layout()
plt.grid(True)
plt.show()

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import numpy as np

# Measured values (in F)
connection_measurement = {
    'C_ABvsCDE': 1.257189992699E-10,
    'C_ABCDvsE': 3.581451055984E-11,
    'C_ABEvsCd': 1.036613057794E-10,
    'C_AvsBCDE': 3.515076277291E-11,
    'C_BvsACDE': 7.638724988931E-11,
    'C_CvsABDE': 7.249036736151E-11,
    'C_DvsABCE': 2.855875033066E-11,
    'C_ACvsBDE': 1.220712831988E-10,
    'C_ADvsBCE': 8.111945395629E-11,
    'C_BC_ADE': 6.237211347932E-11,
}

# Define connection sums (simulation)
connection_sums = {
    'C_ABvsCDE': lambda C: C[2] + C[3] + C[4] + C[5] + C[6] + C[7],
    'C_ABCDvsE': lambda C: C[6] + C[7] + C[8] + C[9],
    'C_ABEvsCd': lambda C: C[2] + C[3] + C[4] + C[5] + C[8] + C[9],
    'C_AvsBCDE': lambda C: C[0] + C[3] + C[5] + C[6],
    'C_BvsACDE': lambda C: C[0] + C[2] + C[4] + C[7],
    'C_CvsABDE': lambda C: C[1] + C[3] + C[4] + C[8],
    'C_DvsABCE': lambda C: C[1] + C[2] + C[5] + C[9],
    'C_ACvsBDE': lambda C: C[0] + C[1] + C[4] + C[5] + C[6] + C[8],
    'C_ADvsBCE': lambda C: C[0] + C[1] + C[2] + C[3] + C[6] + C[9],
    'C_BC_ADE': lambda C: C[0] + C[1] + C[2] + C[3] + C[7] + C[8],
}

# Simulated energy data (We_case_1)
We_case_1 = np.array([
    2.430180364410492e-11, 3.963870021343636e-11, 4.628271783338096e-11,
    1.705897537969348e-11, 5.084215412258796e-11, 5.061364115385603e-11,
    6.041130951808294e-11, 1.521894917175105e-10, 5.994497323888096e-11,
    7.127667730542149e-11
]).reshape((10, 1))

M_case_1 = np.array([
    [1, 0, 0, 1, 0, 1, 1, 0, 0, 0],
    [0, 1, 0, -1, 1, 0, 0, 0, 1, 0],
    [0, 0, 1, -1, 1, -1, 0, 0, 1, 1],
    [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
    [1, 1, 0, 0, 1, 1, 1, 0, 1, 0],
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
    [1, 0, 0, 1, 0, 1, 2, 1, 1, 1],
    [0, 1, 1, -2, 2, -1, 0, 0, 2, 1],
    [0, 1, 0, -1, 1, 0, 1, 1, 2, 1],
    [0, 0, 1, -1, -1, 1, 1, 1, 2, 2],
])

# Compute capacitances
C_case1 = np.linalg.inv(M_case_1 ** 2).dot(2 * We_case_1).flatten()

# Compute simulated and measured capacitance sums in pF
calculated_sums_pf = [connection_sums[key](C_case1) * 1e12 for key in connection_sums]
measured_sums_pf = [val * 1e12 for val in connection_measurement.values()]

# Clean x-axis labels
connection_labels_clean = [
    "AB vs CDE", "ABCD vs E", "ABE vs CD", "A vs BCDE", "B vs ACDE",
    "C vs ABDE", "D vs ABCE", "AC vs BDE", "AD vs BCE", "BC vs ADE"
]

# Plot
plt.figure(figsize=(14, 7))
x_indices = np.arange(len(connection_labels_clean))

plt.scatter(x_indices, calculated_sums_pf, label="simulated capacitance", color='blue', marker='o')
plt.scatter(x_indices, measured_sums_pf, label="Measured capacitance", color='red', marker='x')

# Annotate error ratios
for i, (calc, meas) in enumerate(zip(calculated_sums_pf, measured_sums_pf)):
    ratio = calc / meas
    plt.text(i, calc, f"{ratio:.2f}×", fontsize=9, ha='center', va='bottom', color='blue')

# Formatting
plt.xticks(x_indices, connection_labels_clean, rotation=45)
plt.ylabel('Capacitance (pF)')
plt.title('Comparison of simulated and measured capacitance')
plt.legend(loc='upper right')
plt.tight_layout()
plt.grid(True)

# Save as PDF
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\inductor_42\transformer_58_50_error.pdf"
plt.savefig(pdf_path)


stored_energy_AvsB =  4.40599880490035E-08
C_AvsB_float = (2 * stored_energy_AvsB) /( 57**2)
print(C_AvsB_float)
stored_energy_CvsD = 3.28323539241112E-08
C_CvsD_float = 2 * stored_energy_CvsD / 49**2
print(C_CvsD_float)
C_Avs_B_Model = C_case1[0] + ((C_case1[6] * C_case1[7]) / (C_case1[6] + C_case1[7])) + ((C_case1[3] * C_case1[4]) / (C_case1[3] + C_case1[4])) + ((C_case1[5] * C_case1[2]) / (C_case1[5] + C_case1[2])) + ((C_case1[1] * C_case1[8] * C_case1[9]) / (C_case1[1]  + C_case1[8] + C_case1[9]))
C_Bvs_C_Model = C_case1[1] + ((C_case1[8] * C_case1[9]) / (C_case1[8] + C_case1[9])) + ((C_case1[3] * C_case1[4]) / (C_case1[3] + C_case1[4])) + ((C_case1[5] * C_case1[2]) / (C_case1[5] + C_case1[2])) + ((C_case1[0] * C_case1[8] * C_case1[9]) / (C_case1[0]  + C_case1[8] + C_case1[9]))
C_AvsB_try = C_case1[0]  + ((C_case1[6] * C_case1[7]) / (C_case1[6] + C_case1[7])) + C_case1[2] + C_case1[3] + C_case1[4] + C_case1[5] + C_case1[6] + C_case1[1] + ((C_case1[8] * C_case1[9]) / (C_case1[8] + C_case1[9]))
print(C_Avs_B_Model)
print(C_Bvs_C_Model)
print(C_AvsB_try)

n_sym = 58 / 50
C1 = C_case1[0]
C2 = C_case1[1]
C3 = C_case1[2]
C4 = C_case1[3]
C5 = C_case1[4]
C6 = C_case1[5]
C7 = C_case1[6]
C8 = C_case1[7]
C9 = C_case1[8]
C10 = C_case1[9]

connection_measurement = {
    'C_AvsB (CD open)': 5.636524693E-11,
    'C_AvsB (CD short)': 2.423008688E-11
}



# Calculate C_Meas
C_Meas_open = (
    ((C2*C3*C7 + C2*C3*C8 + C2*C4*C7 + C2*C3*C9 + C2*C4*C8 + C2*C5*C7 + C3*C4*C7 + C2*C3*C10 + C2*C4*C9 + C2*C5*C8 + C2*C6*C7 + C3*C4*C8 + C3*C5*C7 + C2*C4*C10
      + C2*C5*C9 + C2*C6*C8 + C3*C4*C9 + C3*C5*C8 + C2*C5*C10 + C2*C6*C9 + C3*C4*C10 + C3*C5*C9 + C4*C6*C7 + C2*C6*C10 + C2*C7*C9 + C3*C5*C10 + C4*C6*C8 + C5*C6*C7
      + C2*C7*C10 + C2*C8*C9 + C3*C7*C9 + C4*C6*C9 + C5*C6*C8 + C2*C8*C10 + C3*C8*C9 + C4*C6*C10 + C5*C6*C9 + C4*C7*C10 + C5*C6*C10 + C3*C9*C10 + C4*C8*C10 + C5*C7*C10
      + C6*C7*C9 + C4*C9*C10 + C5*C8*C10 + C6*C8*C9 + C5*C9*C10 + C6*C9*C10 + C7*C9*C10 + C8*C9*C10)*n_sym**2)/(C3*C7 + C3*C8 + C4*C7 + C3*C9 + C4*C8 + C5*C7 + C3*C10 + C4*C9
         + C5*C8 + C6*C7 + C4*C10 + C5*C9 + C6*C8 + C5*C10 + C6*C9 + C6*C10 + C7*C9 + C7*C10 + C8*C9 + C8*C10) - ((2*C3*C4*C7 + 2*C3*C4*C8 + 2*C3*C4*C9 + 2*C3*C4*C10
         - 2*C5*C6*C7 + 2*C3*C7*C9 - 2*C5*C6*C8 - 2*C5*C6*C9 - 2*C5*C6*C10 + 2*C4*C8*C10 - 2*C5*C7*C10 - 2*C6*C8*C9)*n_sym)/(C3*C7 + C3*C8 + C4*C7 + C3*C9 + C4*C8 + C5*C7
           + C3*C10 + C4*C9 + C5*C8 + C6*C7 + C4*C10 + C5*C9 + C6*C8 + C5*C10 + C6*C9 + C6*C10 + C7*C9 + C7*C10 + C8*C9 + C8*C10) + (C1*C3*C7 + C1*C3*C8 + C1*C4*C7 + C1*C3*C9 + C1*C4*C8
            + C1*C5*C7 + C1*C3*C10 + C1*C4*C9 + C1*C5*C8 + C1*C6*C7 + C3*C4*C7 + C1*C4*C10 + C1*C5*C9 + C1*C6*C8 + C3*C4*C8 + C1*C5*C10 + C1*C6*C9 + C3*C4*C9 + C3*C6*C7 + C4*C5*C7 + C1*C6*C10
            + C1*C7*C9 + C3*C4*C10 + C3*C6*C8 + C4*C5*C8 + C1*C7*C10 + C1*C8*C9 + C3*C6*C9 + C3*C7*C8 + C4*C5*C9 + C5*C6*C7 + C1*C8*C10 + C3*C6*C10 + C3*C7*C9 + C4*C5*C10 + C4*C7*C8 + C5*C6*C8
            + C3*C7*C10 + C5*C6*C9 + C5*C7*C8 + C4*C8*C9 + C5*C6*C10 + C5*C7*C9 + C6*C7*C8 + C4*C8*C10 + C5*C7*C10 + C6*C8*C9 + C6*C8*C10 + C7*C8*C9 + C7*C8*C10)/(C3*C7 + C3*C8 + C4*C7 + C3*C9
             + C4*C8 + C5*C7 + C3*C10 + C4*C9 + C5*C8 + C6*C7 + C4*C10 + C5*C9 + C6*C8 + C5*C10 + C6*C9 + C6*C10 + C7*C9 + C7*C10 + C8*C9 + C8*C10)
)
print(C_Meas_open)

C_Meas_short = (
    C1 + (C3*C4*C7 + C3*C4*C8 + C3*C4*C9 + C3*C6*C7 + C4*C5*C7 + C3*C4*C10 + C3*C6*C8 + C4*C5*C8 + C3*C6*C9 + C3*C7*C8 + C4*C5*C9 + C5*C6*C7 + C3*C6*C10 + C3*C7*C9
          + C4*C5*C10 + C4*C7*C8 + C5*C6*C8 + C3*C7*C10 + C5*C6*C9 + C5*C7*C8 + C4*C8*C9 + C5*C6*C10 + C5*C7*C9 + C6*C7*C8 + C4*C8*C10 + C5*C7*C10 + C6*C8*C9
          + C6*C8*C10 + C7*C8*C9 + C7*C8*C10)/(C3*C7 + C3*C8 + C4*C7 + C3*C9 + C4*C8 + C5*C7 + C3*C10 + C4*C9 + C5*C8 + C6*C7 + C4*C10 + C5*C9 + C6*C8 + C5*C10
                                               + C6*C9 + C6*C10 + C7*C9 + C7*C10 + C8*C9 + C8*C10)
)
print(C_Meas_short)


den = C3*C7 + C3*C8 + C4*C7 + C3*C9 + C4*C8 + C5*C7 + C3*C10 + C4*C9 + \
      C5*C8 + C6*C7 + C4*C10 + C5*C9 + C6*C8 + C5*C10 + C6*C9 + C6*C10 + \
      C7*C9 + C7*C10 + C8*C9 + C8*C10

num1 = (C2*C3*C7 + C2*C3*C8 + C2*C4*C7 + C2*C3*C9 + C2*C4*C8 + C2*C5*C7 +
        C3*C4*C7 + C2*C3*C10 + C2*C4*C9 + C2*C5*C8 + C2*C6*C7 + C3*C4*C8 +
        C3*C5*C7 + C2*C4*C10 + C2*C5*C9 + C2*C6*C8 + C3*C4*C9 + C3*C5*C8 +
        C2*C5*C10 + C2*C6*C9 + C3*C4*C10 + C3*C5*C9 + C4*C6*C7 + C2*C6*C10 +
        C2*C7*C9 + C3*C5*C10 + C4*C6*C8 + C5*C6*C7 + C2*C7*C10 + C2*C8*C9 +
        C3*C7*C9 + C4*C6*C9 + C5*C6*C8 + C2*C8*C10 + C3*C8*C9 + C4*C6*C10 +
        C5*C6*C9 + C4*C7*C10 + C5*C6*C10 + C3*C9*C10 + C4*C8*C10 + C5*C7*C10 +
        C6*C7*C9 + C4*C9*C10 + C5*C8*C10 + C6*C8*C9 + C5*C9*C10 + C6*C9*C10 +
        C7*C9*C10 + C8*C9*C10) * n_sym**2

num2 = (-2*C3*C4*C7 - 2*C3*C4*C8 - 2*C3*C4*C9 - 2*C3*C4*C10 + 2*C5*C6*C7 +
        2*C3*C7*C9 - 2*C5*C6*C8 - 2*C5*C6*C9 - 2*C5*C6*C10 + 2*C4*C8*C10 -
        2*C5*C7*C10 - 2*C6*C8*C9) * n_sym

num3 = (C1*C3*C7 + C1*C3*C8 + C1*C4*C7 + C1*C3*C9 + C1*C4*C8 + C1*C5*C7 +
        C1*C3*C10 + C1*C4*C9 + C1*C5*C8 + C1*C6*C7 + C3*C4*C7 + C1*C4*C10 +
        C1*C5*C9 + C1*C6*C8 + C3*C4*C8 + C1*C5*C10 + C1*C6*C9 + C3*C4*C9 +
        C3*C6*C7 + C4*C5*C7 + C1*C6*C10 + C1*C7*C9 + C3*C4*C10 + C3*C6*C8 +
        C4*C5*C8 + C1*C7*C10 + C1*C8*C9 + C3*C6*C9 + C3*C7*C8 + C4*C5*C9 +
        C5*C6*C7 + C1*C8*C10 + C3*C6*C10 + C3*C7*C9 + C4*C5*C10 + C4*C7*C8 +
        C5*C6*C8 + C3*C7*C10 + C5*C6*C9 + C5*C7*C8 + C4*C8*C9 + C5*C6*C10 +
        C5*C7*C9 + C6*C7*C8 + C4*C8*C10 + C5*C7*C10 + C6*C8*C9 + C6*C8*C10 +
        C7*C8*C9 + C7*C8*C10)

C_meas_open_sym = (num1 - num2 + num3) / den
print(C_meas_open_sym)

import numpy as np
import matplotlib.pyplot as plt

# Define symbolic or previously calculated expressions
# These are assumed to be evaluated earlier in your script:
# C_Meas_open and C_Meas_short

# Use symbolic results (or numerically evaluated ones)
import numpy as np
import matplotlib.pyplot as plt

# Simulated values calculated from equations
simulated_vals = [C_Meas_open, C_Meas_short]

# Measured values in Farads
measured_vals = [5.636524693E-11, 2.423008688E-11]

# Labels and positions
labels = ["A vs B\n(CD open)", "A vs B\n(CD short)"]  # Add line break to reduce width
x = np.arange(len(labels))

# Plot
plt.figure(figsize=(6, 5))  # Smaller width for tighter space
plt.scatter(x, simulated_vals, label="Simulated", marker='o', color='blue')
plt.scatter(x, measured_vals, label="Measured", marker='x', color='red')

# Annotate error ratio
for i, (sim, meas) in enumerate(zip(simulated_vals, measured_vals)):
    ratio = sim / meas
    plt.text(i, sim, f"{ratio:.2f}×", ha='center', va='bottom', color='blue')

# Customize layout
plt.xticks(x, labels)
plt.ylabel("Capacitance (F)")
plt.title("Comparison of Simulated and Measured Capacitance\nA vs B (CD open/short)", fontsize=12)
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save the plot
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\inductor_42\transformer_open_short_comparison_compact.pdf"
plt.savefig(pdf_path)
plt.show()

# Simulated and measured values
simulated_vals = [C_Meas_open, C_Meas_short]
measured_vals = [5.636524693E-11, 2.423008688E-11]

# Labels
labels = ["A vs B (CD open)", "A vs B (CD short)"]
y_pos = np.arange(len(labels))

# Figure setup
plt.figure(figsize=(12, 4))
bar_height = 0.3

# Bars
plt.barh(y_pos - bar_height/2, simulated_vals, height=bar_height, color='tab:blue', label="Simulated")
plt.barh(y_pos + bar_height/2, measured_vals, height=bar_height, color='tab:red', label="Measured")

# Annotate error ratios
for i, (sim, meas) in enumerate(zip(simulated_vals, measured_vals)):
    ratio = sim / meas
    plt.text(sim * 1.01, i - bar_height/2, f"{ratio:.2f}×", va='center', fontsize=9, color='blue')

# Layout and labels
plt.yticks(y_pos, labels)
plt.xlabel("Capacitance (F)")
plt.title("Comparison of simulated and measured capacitance\nA vs B (CD open/short)", fontsize=10)
plt.grid(axis='x', linestyle='--', alpha=0.6)
plt.legend(loc='upper right')
plt.tight_layout()

# Save
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\inductor_42\transformer_open_short_bar_centered.pdf"
plt.savefig(pdf_path)
plt.show()


# Convert values to picoFarads
simulated_vals_pf = [C_Meas_open * 1e12, C_Meas_short * 1e12]
measured_vals_pf = [5.636524693E-11 * 1e12, 2.423008688E-11 * 1e12]

# Labels
labels = ["A vs B (CD open)", "A vs B (CD short)"]
y_pos = np.arange(len(labels))

# Figure setup
plt.figure(figsize=(12, 4))
bar_height = 0.3

# Bars
plt.barh(y_pos - bar_height/2, simulated_vals_pf, height=bar_height, color='tab:blue', label="Simulated")
plt.barh(y_pos + bar_height/2, measured_vals_pf, height=bar_height, color='tab:red', label="Measured")

# Annotate error ratios
for i, (sim, meas) in enumerate(zip(simulated_vals_pf, measured_vals_pf)):
    ratio = sim / meas
    plt.text(sim * 1.01, i - bar_height/2, f"{ratio:.2f}×", va='center', fontsize=9, color='blue')

# Layout and labels
plt.yticks(y_pos, labels)
plt.xlabel("Capacitance (pF)")
plt.title("Comparison of simulated and measured capacitance\nA vs B (CD open/short)", fontsize=10)
plt.grid(axis='x', linestyle='--', alpha=0.6)

# Legend
# plt.legend(loc='upper right', bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
plt.legend(loc='upper right')

plt.tight_layout()

# Save
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\inductor_42\transformer_open_short_bar_centered_pf.pdf"
plt.savefig(pdf_path)
plt.show()
