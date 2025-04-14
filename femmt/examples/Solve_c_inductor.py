import math

import numpy as np
import matplotlib.pyplot as plt
from sympy.codegen import Print

# Define the M matrix based on your corrected scenarios

M_case_1 = np.array([
    [1, 1, 0],  # Scenario 1
    [0, 1, 1],  # Scenario 2
    #[1, 0.4327, 0.5673],  # Scenario 3
    [1, 0, 1]
    # [1, 2 ,1],
    # [1, 0.5, -0.5]
])
# M_case_1 = np.array([
#     [1, 0],  # Scenario 1
#     [0, 1],  # Scenario 2
#     #[1, 0.4327, 0.5673],  # Scenario 3
#     # [0, 1 ,1]
# ])

M_squared1 = M_case_1 ** 2

# We_case_1 = np.array([
#     3.25739581637384E-12,  # Scenario 1
#     2.28651876671513E-12,   # Scenario 2
#     #1.467349670661194e-12,
#     7.330786489864692e-12# Scenario 3
# ])
We_case_1 = np.array([
    5.23950306321443e-12,  # Scenario 1
    1.420017752187458e-11, # Scenario 2
    5.05540791441735e-12,# Scenario 3
])
We_case1 = We_case_1.reshape((3, 1))
det_M1 = np.linalg.det(M_squared1)
print(f"Determinant of M: {det_M1}")

if det_M1 == 0:
    print("Matrix M is singular and cannot be inverted. Consider revising your simulation scenarios.")
else:
    C_case1 = 2 * np.linalg.inv(M_squared1).dot(We_case1)

    # Flatten C to a 1D array for easier indexing
    C_case1 = C_case1.flatten()

    # Display
    capacitance_pairs = ["A-B", "A-0", "B-0", "A-C", "B-C",
                         "A-D", "A-E", "B-E", "C-E", "D-E"]
    print("\nIndividual Capacitances:")
    for idx, capacitance in enumerate(C_case1, start=1):
        print(f"C{idx} ({capacitance_pairs[idx-1]}): {capacitance:.5e} F")

    A = 1.594250846194489e-12 * 2
    print(A)
    B = ((A-2.80E-12)/A)*100
    print(B)
    # C = 1.43843e-11 - 3.90527e-12
    # print(C)

    Ctt = 2.213628432187502e-12
    Cts = 2.166057416323814e-12

    C_ab = Ctt + Cts / 2  # Initial for n = 2

    # for n in range(4, 44, 2):  # 4 to 42
    #     C_ab = (C_ab * (Ctt / 2)) / (C_ab + (Ctt / 2)) + (Cts / 2)
    c_AB = 2.213628432187502e-12 * 1.336
    print(f"C_AB(42) = {C_ab:.5e} F")

    C_AB_Till = -3.91159e-12 + (1.43906e-11 * 1.40098e-11 / (1.43906e-11 + 1.40098e-11))
    print("C_Till:", C_AB_Till)

    C_AB_Othman = (-3.91159e-12 + 0.25 * 1.43906e-11 + 0.25 * 1.40098e-11)
    print("C_Othman:", C_AB_Othman)

    print("A vs BE:", -3.91159e-12 + 1.43906e-11)
    print("B vs AE:", -3.91159e-12 + 1.40098e-11)
    print("AB vs E:", 1.43906e-11 + 1.40098e-11)
    # m = 1.094475350790898e-12 * 1.336
    # print(m)
    #
    # n = 1.162248224277244e-09 * 2 / (14*14)
    # print(n)
    # # inductor 81
    # C_AB_20 = 2 * 1.022199231318327e-12
    # C_tt_41 = 2.096884260266433e-12
    # C_AB_42 = 2 * 1.594250846194489e-12
    # C_AB_62 = 4.533100005535176e-12 * 2
    # C_AB_81 = 1.838169938418392e-11 * 2
    # C_AB_102 = 1.379046892383806e-11 * 2
    # C_AB_168 = 1.523104163816231e-11 * 2
    # C_tt_200 = 1.751260449573454e-12
    # C_AB_200 = 1.282991822772131e-11 * 2
    #
    # # test
    # d = 1e-3
    # C_tt_5 = 6.811551995306939e-13
    # C_AB_5 = 7.7571660608902e-13 * 2
    # C_tt_7 = 6.811345919870133e-13 * 1.336
    # C_AB_7 = 8.436371098605758e-13 * 2
    # C_tt_14 = 6.811324381178942e-13 * 1.336
    # C_AB_14 = 1.168419425373008e-12 *2
    # C_tt_17 = 6.811324603992539e-13 * 1.336
    # C_AB_17 = 1.341620619464279e-12 * 2
    # #
    # d = 5e-4
    # C_tt_7 = 1.068186865777103e-12 * 1.336
    # C_AB_5 = 7.933065628771767e-13 * 2
    # C_tt_7 = 1.068186865777103e-12 * 1.336
    # C_AB_7 = 8.144139085584939e-13 * 2
    # C_tt_14 = 1.068176538193686e-12 * 1.336
    # C_AB_14 = 1.002974472605233e-12 * 2
    # C_tt_17 = 1.068176592940226e-12 * 1.336
    # C_Ab_17 = 1.11459153041314e-12 * 2
    # #
    # d = 2.5e-4
    # C_tt_5 = 1.522423785044244e-12 * 1.336
    # C_AB_5 = 8.416078977508367e-13 * 2
    # C_tt_7 = 1.522260009119674e-12 * 1.336
    # C_AB_7 = 8.242435453777094e-13 * 2
    # C_AB_14 = 9.433097864360272e-13 * 2
    # C_tt_14 = 1.522260009119674e-12 * 1.336
    # C_Ab_17 = 1.00777334598581e-12 * 2
    # C_tt_17 = 1.522260009119674e-12 * 1.336
    # #
    # d = 1.25e-4
    # A_AB_5 = 9.120132866389397e-13 * 2
    # C_AB_7 = 8.604810238430487e-13 * 2
    # C_AB_14 = 9.181645383291066e-13 * 2
    # C_AB_17 = 9.729401889377677e-13 * 2
    # A_AB_20 = 1.026168261308761e-12 * 2
    # C_AB_30 = 4.40483641299e-13 * 2
    # C_tt_40 = 2.246371038719414e-12 * 1.336
    # C_AB_40 = 1.53383569696869e-12 * 2
    #
    # d = 6.25e-6
    # C_tt_50 = 8.393292613868711e-12 * 1.336
    # C_tt_40 = 8.393256567924022e-12 * 1.336
    # C_tt_30 = 8.393256581649125e-12 * 1.336
    # C_tt_70 = 8.393256581649125e-12 * 1.336
    # C_AB_55 = 1.680063320492602e-12 * 2
    # C_Ab_40 = 1.360836592462348e-12 * 2
    # C_AB_30 = 1.195789728976239e-12 * 2
    # C_AB_7 = 1.388464814610847e-12 * 2
    # C_AB_5 = 1.718302418283709e-12 * 2

    # M = 0.95  * math.sqrt(0.001531946165 + 0.001141382479)
    # print(M)
    # L_1 = 0.001531946165
    # L_2 = 0.001141382479
    # K = 0.9972568281187678
    #
    # # Calculate Mutual Inductance M
    # M = K * math.sqrt(L_1 * L_2)
    #
    # # Calculate Lm
    # L_m = L_1 - (M ** 2 / L_2)
    #
    # print(L_m)





