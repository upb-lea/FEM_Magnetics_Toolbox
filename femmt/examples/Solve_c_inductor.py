
import numpy as np

# Define the M matrix based on your corrected scenarios

M_case_1 = np.array([
    [1, 1, 0],  # Scenario 1
    [0, 1, 1],  # Scenario 2
    [1, 0, 1]  # Scenario 3
])

M_squared1 = M_case_1 ** 2
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

    # A = 1.594250846194489e-12 * 2
    # print(A)
    # B = ((A-2.80E-12)/A)*100
    # print(B)
    # # C = 1.43843e-11 - 3.90527e-12
    # # print(C)
    #
    # Ctt = 2.213628432187502e-12
    # Cts = 2.166057416323814e-12
    #
    # C_ab = Ctt + Cts / 2  # Initial for n = 2
    #
    # # for n in range(4, 44, 2):  # 4 to 42
    # #     C_ab = (C_ab * (Ctt / 2)) / (C_ab + (Ctt / 2)) + (Cts / 2)
    # c_AB = 2.213628432187502e-12 * 1.336
    # print(f"C_AB(42) = {C_ab:.5e} F")
    #
    # C_AB_Till = -3.91159e-12 + (1.43906e-11 * 1.40098e-11 / (1.43906e-11 + 1.40098e-11))
    # print("C_Till:", C_AB_Till)
    #
    # C_AB_Othman = (-3.91159e-12 + 0.25 * 1.43906e-11 + 0.25 * 1.40098e-11)
    # print("C_Othman:", C_AB_Othman)
    #
    # print("A vs BE:", -3.91159e-12 + 1.43906e-11)
    # print("B vs AE:", -3.91159e-12 + 1.40098e-11)
    # print("AB vs E:", 1.43906e-11 + 1.40098e-11)




