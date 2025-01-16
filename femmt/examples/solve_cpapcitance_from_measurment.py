import numpy as np

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
# Define the energy vector (replace with your actual measured energies)

# with error 50%
# adjusting
# We_no_bobbin = np.array([
#     1.107533022806718e-11,  # Scenario 1
#     2.998536727832599e-11,   # Scenario 2
#     2.441518883756004e-11,   # Scenario 3
#     1.126487136125143e-11,  # Scenario 4
#     3.436564097908541e-11,  # Scenario 5
#     3.11986239282366e-11,  # Scenario 6
#     2.014172856588652e-11,  # Scenario 7
#     1.534054057620717e-11,  # Scenario 8
#     2.924182091720404e-11,  # Scenario 9
#     3.222388455288184e-11,   # Scenario 10
#     # 1.114304711308078e-11
# ])
# We_with_bobbin = np.array([
#     2.358792372177351e-11,  # Scenario 1
#     3.10183406267461e-11,   # Scenario 2
#     2.56393791428112e-11,   # Scenario 3
#     1.248543827275428e-11,  # Scenario 4
#     5.172790534894289e-11,  # Scenario 5
#     4.484317223246496e-11,  # Scenario 6
#     3.387134997216844e-11,  # Scenario 7
#     1.752841658974935e-11,  # Scenario 8
#     3.15035211362564e-11,  # Scenario 9
#     3.478179986637866e-11,   # Scenario 10
#     # 2.350977482015328e-11    # 11
# ])
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

# Ensure We is a column vector
# We = We_no_bobbin.reshape((10, 1))
# We = We_with_bobbin.reshape((10, 1))
We = We_with_modified_bobbin.reshape((10, 1))
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

energy_with_bobbin = 7.676930117930947e-11
energy_without_bobbin = 1.107533022806718e-11

difference = energy_with_bobbin - energy_without_bobbin
percentage_difference = (difference / energy_without_bobbin) * 100

# Print the results
print(f"Absolute difference in energy: {difference:e} joules")
print(f"Percentage difference: {percentage_difference:.2f}%")
