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

# We_case_1 = np.array([
#     1.066244380397435e-11,  # Scenario 1
#     1.193107829042975e-11,   # Scenario 2
#     3.008316244330874e-11,   # Scenario 3
#     1.114015694663501e-11,  # Scenario 4
#     # 5.805224198514209e-12,  # Scenario 5
#     1.489056178747461e-11,
#     1.72746715995602e-11,  # Scenario 6
#     2.714628935917936e-11,  # Scenario 7
#     7.311772673214487e-11,  # Scenario 8
#     # 3.018397739929745e-11,  # Scenario 9
#     2.8000140286046e-11,
#     5.428529015625401e-11,   # Scenario 10
#     # 3.044751388584012e-11    # 11
# ])
# increase distance between primary and sec
# We_case_1 = np.array([
#     9.596271818853584e-12,  # Scenario 1
#     9.714731614410892e-12,   # Scenario 2
#     2.344374039374051e-11,   # Scenario 3
#     2.219089948626731e-11,  # Scenario 4
#     # 5.805224198514209e-12,  # Scenario 5
#     1.456000062094578e-11,
#     1.974805257372915e-11,  # Scenario 6
#     4.437799160359449e-11,  # Scenario 7
#     5.74256981261533e-11,  # Scenario 8
#     # 3.018397739929745e-11,  # Scenario 9
#     4.063811656324148e-11,
#     6.58710675473648e-11,   # Scenario 10
#     # 3.044751388584012e-11    # 11
# ])
# winding direction is different
We_case_1 = np.array([
    1.319424236783797e-11,  # Scenario 1
    1.162485250286133e-11,   # Scenario 2
    2.956127865036595e-11,   # Scenario 3
    1.681958299078067e-11,  # Scenario 4
    # 5.805224198514209e-12,  # Scenario 5
    7.551174869346277e-12,
    1.858518120804068e-11,  # Scenario 6
    4.212456498984411e-11,  # Scenario 7
    7.162414623050972e-11,  # Scenario 8
    # 3.018397739929745e-11,  # Scenario 9
    3.417719489857039e-11,
    5.685486200502775e-11,   # Scenario 10
    # 3.044751388584012e-11    # 11
])



We_case1 = We_case_1.reshape((10, 1))
det_M1 = np.linalg.det(M_squared1)
print(f"Determinant of M: {det_M1}")

if det_M1 == 0:
    print("Matrix M is singular and cannot be inverted. Consider revising your simulation scenarios.")
else:
    C_case1 = 2 * np.linalg.inv(M_squared1).dot(We_case1)

    # Flatten C to a 1D array for easier indexing
    C_case1 = C_case1.flatten()

    # Display
    # capacitance_pairs = ["A-B", "C-D", "A-D", "B-C", "B-D",
    #                      "A-C", "A-E", "B-E", "D-E", "C-E"]
    capacitance_pairs = ["A-B", "C-D", "B-C", "A-D", "B-D",
                         "A-C", "A-E", "B-E", "D-E", "C-E"]
    print("\nIndividual Capacitances:")
    for idx, capacitance in enumerate(C_case1, start=1):
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

