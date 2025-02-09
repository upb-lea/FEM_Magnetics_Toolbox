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
We_with_modified_ins_cond_epsilon_r_is_one = np.array([
    1.896646196533772e-11,  # Scenario 1
    1.785239296253184e-11,   # Scenario 2
    1.818042720391898e-11,   # Scenario 3
    1.703657012251492e-11,  # Scenario 4
    5.295722229338796e-11,  # Scenario 5
    3.432119340610629e-12,  # Scenario 6
    2.079145395412682e-11,  # Scenario 7
    2.104150832013626e-11,  # Scenario 8
    3.117737033725389e-12,  # Scenario 9
    5.070835812665073e-11,   # Scenario 10
    #3.044751388584012e-11    # 11
])
# femm
We_with_modified_ins_cond_epsilon_r_is_one_femm = np.array([
    1.904312842747512e-11,  # Scenario 1
    1.792092892424583e-11,   # Scenario 2
    1.825802808810904e-11,   # Scenario 3
    1.710347501700135e-11,  # Scenario 4
    5.316169679438457e-11,  # Scenario 5
    3.438355010228581e-12,  # Scenario 6
    2.087385227086202e-11,  # Scenario 7
    2.11273133020259e-11,  # Scenario 8
    3.120147640149721e-12,  # Scenario 9
    5.091331977408417e-11,   # Scenario 10
    #3.044751388584012e-11    # 11
])
# # with_modified_new_scheme_cond
# We_with_modified_new_scheme_cond = np.array([
#     2.773350214336289e-11,  # Scenario 1
#     2.34209643382356e-11,   # Scenario 2
#     2.679920014773644e-11,   # Scenario 3
#     2.271721133209746e-11,  # Scenario 4
#     7.431795139552544e-11,  # Scenario 5
#     3.518148742480362e-12,  # Scenario 6
#     2.805293125056189e-11,  # Scenario 7
#     2.839120815089925e-11,  # Scenario 8
#     2.986424441328575e-12,  # Scenario 9
#     7.19653566182524e-11,   # Scenario 10
#     #3.044751388584012e-11    # 11
# ])
We_with_modified_ins_cond_epsilon_r_is_three= np.array([
    2.756508863415358e-11,  # Scenario 1
    2.456509071896705e-11,   # Scenario 2
    2.671265230458517e-11,   # Scenario 3
    2.376064534573603e-11,  # Scenario 4
    7.501676716808256e-11,  # Scenario 5
    3.527239633582984e-12,  # Scenario 6
    2.929880303521695e-11,  # Scenario 7
    2.957319959896681e-11,  # Scenario 8
    3.185381956923266e-12,  # Scenario 9
    7.272732069923849e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])
We_with_correct_distances = np.array([
    1.403087465503489e-11,  # Scenario 1
    1.357402345230772e-11,   # Scenario 2
    1.10915378827943e-11,   # Scenario 3
    1.052400658622305e-11,  # Scenario 4
    3.966934784395457e-11,  # Scenario 5
    5.529572561534309e-12,  # Scenario 6
    1.558611820075e-11,  # Scenario 7
    1.582570545902131e-11,  # Scenario 8
    5.174954071864558e-12,  # Scenario 9
    3.083036233763695e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])
We_with_correct_distances_with_right_ins = np.array([
    1.931909816212589e-11,  # Scenario 1
    1.830944905932953e-11,  # Scenario 2
    1.634784967968647e-11,  # Scenario 3
    1.52939928649424e-11,  # Scenario 4
    5.43011404953074e-11,  # Scenario 5
    5.511522514843165e-12,  # Scenario 6
    2.099970806544051e-11,  # Scenario 7
    2.12423267742129e-11,  # Scenario 8
    5.163966193659957e-12,  # Scenario 9
    4.548916618682985e-11,  # Scenario 10
    # 3.044751388584012e-11    # 11
])
We_with_correct_distances_with_right_ins_with_epsilon_r_of_bobbin_four = np.array([
    1.876090034312706e-11,  # Scenario 1
    1.781887585976625e-11,   # Scenario 2
    1.618423702929472e-11,   # Scenario 3
    1.525025261528419e-11,  # Scenario 4
    5.291013528769602e-11,  # Scenario 5
    4.898009590389586e-12,  # Scenario 6
    2.039885523959554e-11,  # Scenario 7
    2.058946412878626e-11,  # Scenario 8
    4.640690923584991e-12,  # Scenario 9
    4.52774739951255e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])
We_Till = np.array([
    2.018229492973006e-11,  # Scenario 1
    1.719891768514534e-11,   # Scenario 2
    1.620818447105329e-11,   # Scenario 3
    1.558863170201864e-11,  # Scenario 4
    5.521783962238933e-11,  # Scenario 5
    5.679027828644863e-12,  # Scenario 6
    1.995202333129267e-11,  # Scenario 7
    2.056192892148821e-11,  # Scenario 8
    4.911796464153383e-12,  # Scenario 9
    4.57584432259361e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])
# # Ensure We is a column vector
# We = We_no_bobbin.reshape((10, 1))
# We = We_with_bobbin.reshape((10, 1))
# We = We_with_modified_bobbin.reshape((10, 1))
We = We_Till.reshape((10, 1))
# Check if M is invertible by computing its determinant
det_M = np.linalg.det(M_squared)
print(f"Determinant of M: {det_M}")

if det_M == 0:
    print("Matrix M is singular and cannot be inverted. Consider revising your simulation scenarios.")
else:
    # Solve for C: C = M^{-1} * (2 * We)
    # C = np.linalg.inv(M_squared).dot(2 * We)
    C = 2 * np.linalg.inv(M_squared).dot(We)
    # C = 2 * M_squared.dot(We)

    # Flatten C to a 1D array for easier indexing
    C = C.flatten()

    # Display
    capacitance_pairs = ["A-B", "C-D", "B-D", "A-C", "B-C",
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
energy_data = [We_with_modified_ins_cond_epsilon_r_is_one, We_with_modified_ins_cond_epsilon_r_is_one_femm , We_with_modified_ins_cond_epsilon_r_is_three, We_with_correct_distances, We_with_correct_distances_with_right_ins, We_with_correct_distances_with_right_ins_with_epsilon_r_of_bobbin_four, We_Till]
energy_labels = ["modified_ins_cond_epsilon_r_is_one", "modified_ins_cond_epsilon_r_is_one_femm", "modified_ins_cond_epsilon_r_is_three", "We_with_correct_ditances", "with_correct_distances_with_right_ins", "with_correct_distances_with_right_ins_with_epsilon_r_of_bobbin_four", "Till"]

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
# print

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

condition_number = np.linalg.cond(M_squared)
print(f"Condition number of M_squared: {condition_number}")

# # Define a function to calculate and print connection sums
# def print_and_plot_sums(C, label):
#     print(f"\n--- Connection Sums for Simulation: {label} ---")
#     for name, func in connection_sums.items():
#         sum_value = func(C)
#         measured_value = connection_measurement[name]
#         print(f"{name}: {sum_value:.5e} F (Measured: {measured_value:.5e} F)")
#
# # Iterate over all simulations and calculate sums
# for energy, label in zip(energy_data, energy_labels):
#     We = energy.reshape((10, 1))
#     C = np.linalg.inv(M ** 2).dot(2 * We).flatten()
#     print_and_plot_sums(C, label)
#
#     # Plot the calculated and measured sums
#     calculated_sums = [connection_sums[key](C) for key in connection_sums]
#     measured_sums = list(connection_measurement.values())
#
#     plt.scatter(range(len(calculated_sums)), calculated_sums, label=f"Calculated - {label}", marker='o')
#     plt.scatter(range(len(measured_sums)), measured_sums, label="Measured Connections", color="red", marker='x')
#
# # Set plot details
# plt.xticks(range(len(connection_sums)), list(connection_sums.keys()), rotation=45)
# plt.ylabel('Sum of Capacitances (F)')
# plt.title('Connection Sums for Different Cases and Measured Connections')
# plt.legend()
# plt.grid(True)
# plt.show()
#
#
# Define a function to calculate, print, and plot connection sums with error
def print_and_plot_sums_with_error(C, label):
    print(f"\n--- Connection Sums for Simulation: {label} ---")
    errors = []

    for name, func in connection_sums.items():
        sum_value = func(C)
        measured_value = connection_measurement[name]
        error = abs(sum_value - measured_value) / measured_value * 100
        errors.append(error)

        print(f"{name}: {sum_value:.5e} F (Measured: {measured_value:.5e} F) | Error: {error:.2f}%")

    return errors


# Iterate over all simulations and calculate sums and errors
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

# def Till_procedure():
#     We = [2.018229492973006e-11, 1.719891768514534e-11, 1.620818447105329e-11, 1.558863170201864e-11, 5.521783962238933e-11, 5.679027828644863e-12
#           1.995202333129267e-11, 2.056192892148821e-11, 4.911796464153383e-12, 4.57584432259361e-11]
#     We_1 = 2.018229492973006e-11
#     We_2 = 1.719891768514534e-11
#     We_3 = 1.620818447105329e-11
#     We_4 = 1.558863170201864e-11
#     We_5 = 5.521783962238933e-11
#     We_6 = 5.679027828644863e-12
#     We_7 = 1.995202333129267e-11
#     We_8 = 2.056192892148821e-11
#     We_9 = 4.911796464153383e-12
#     We_10= 4.57584432259361e-11
#
#     C_ABvsCDE_meas = 1.51e-10,
#     C_ABEvsCd_meas = 1.35e-10,
#     C_ABCDvsE_meas = 3.21e-11,
#     C_AvsBCDE_meas = 1.8145e-10,
#     C_BvsACDE_meas = 1.9565e-10,
#     C_CvsABDE_meas= 2.0014e-10,
#     C_DvsABCE_meas= 1.4616e-10,
#     C_ACvsBDE_meas= 1.3501e-11,
#     C_ADvsBCE_meas= 1.3021e-10,
#     C_BC_ADE_meas= 9.0531e-10,
#     C_BD_ACE_meas= 1.3806e-11,
#
#
#     for we in We:
#         C=We(we) * 2 / 1


def Till_procedure():
    We_values = [
        2.018229492973006e-11, 1.719891768514534e-11, 1.620818447105329e-11, 1.558863170201864e-11,
        5.521783962238933e-11, 5.679027828644863e-12, 1.995202333129267e-11, 2.056192892148821e-11,
        4.911796464153383e-12, 4.57584432259361e-11, 2.223643475336065e-11
    ]

    # Measured capacitance
    measured_sums = {
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
        'C_BD_ACE': 1.3806e-11,
    }


    # Simulate the capacitance calculation
    simulated_capacitances = [2 * we / 1 for we in We_values]

    # Compare each simulated value to its respective measured value
    print("\n--- Comparison of Simulated and Measured Capacitance Sums ---")
    for (name, measured_value), simulated_value in zip(measured_sums.items(), simulated_capacitances):
        error = abs(simulated_value - measured_value) / simulated_value * 100
        # error = simulated_value / measured_value * 100
        print(f"{name}: Simulated = {simulated_value:.5e} F, Measured = {measured_value:.5e} F, Error = {error:.2f}%")

    # Plotting the results
    plt.figure(figsize=(12, 6))

    # Plot measured values
    plt.scatter(range(len(measured_sums)), list(measured_sums.values()), color="red", marker='x', label="Measured Capacitance")

    # Plot simulated values
    plt.scatter(range(len(simulated_capacitances)), simulated_capacitances, color="blue", marker='o', label="Simulated Capacitance")

    # Set plot details
    plt.xticks(range(len(measured_sums)), list(measured_sums.keys()), rotation=45)
    plt.ylabel('Capacitance (F)')
    plt.title('Comparison of Simulated and Measured Capacitance Sums')
    plt.legend()
    plt.grid(True)
    plt.show()

    # C_ABvsCDE_new = 2 * (We_values[2] + We_values[3] + We_values[4] + We_values[5] + We_values[6] + We_values[7])
    # print("C_ABvsCDE_new:", C_ABvsCDE_new)
    # # C_ABEvsCD_new = 2 * (We_values[2] + We_values[3] + We_values[4] + We_values[5] + We_values[6] + We_values[7])
    # print("C_ABvsCDE_new:", C_ABvsCDE_new)
    # C_AvsBCDE_new = 2 * (We_values[0] + We_values[3] + We_values[5] + We_values[6])
    # print("C_AvsBCDE_new:", C_AvsBCDE_new)
    # C_BvsACDE_new = 2 * (We_values[0] + We_values[2] + We_values[4] + We_values[7])
    # print("C_BvsACDE_new:", C_BvsACDE_new)
    # C_CvsABDE_new = 2 * (We_values[1] + We_values[3] + We_values[4] + We_values[8])
    # print("C_CvsABDE_new:", C_CvsABDE_new)
    # C_DvsABCE_new = 2 * (We_values[1] + We_values[2] + We_values[5] + We_values[9])
    # print("C_DvsABCE_new:", C_DvsABCE_new)
    # C_BCvsADE_new = 2 * (We_values[0] + We_values[1] + We_values[2] + We_values[3] + We_values[7] + We_values[8])
    # print("C_BCvsADE_new:", C_BCvsADE_new)
    # C_BDvsACE_new = 2 * (We_values[0] + We_values[1] + We_values[2] + We_values[3] + We_values[7] + We_values[8])
    # print("C_BCvsADE_new:", C_BCvsADE_new)


# Call the function to execute the procedure
Till_procedure()










