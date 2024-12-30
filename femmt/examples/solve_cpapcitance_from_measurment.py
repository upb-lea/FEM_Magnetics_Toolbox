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
])



# Define the energy vector (replace with your actual measured energies)
# a little bit ok
# We = np.array([
#     1.06419865439424E-11,  # Scenario 1
#     4.73039528460407E-11,   # Scenario 2
#     2.89205878975124E-11,   # Scenario 3
#     1.40733123145817E-11,  # Scenario 4
#     5.6876634073016E-11,  # Scenario 5
#     3.27093626809716E-11,  # Scenario 6
#     2.03412410765005E-11,  # Scenario 7
#     1.66543080174277E-11,  # Scenario 8
#     2.94952547619938E-11,  # Scenario 9
#     5.45439675817757E-11   # Scenario 10
# ])
# with error 50%
# We = np.array([
#     1.345841649077666e-11,  # Scenario 1
#     3.37605135169375e-11,   # Scenario 2
#     3.043399332689536e-11,   # Scenario 3
#     1.871339220701808e-11,  # Scenario 4
#     3.638210566173711e-11,  # Scenario 5
#     3.893581180112835e-11,  # Scenario 6
#     2.997845514635252e-11,  # Scenario 7
#     2.565655626511367e-11,  # Scenario 8
#     3.640981437360513e-11,  # Scenario 9
#     3.39436495009198e-11   # Scenario 10
# ])
# We = np.array([
#     1.12005406805873E-11,  # Scenario 1
#     2.64859326385962E-11,   # Scenario 2
#     2.41637215441754E-11,   # Scenario 3
#     1.50273531357139E-11,  # Scenario 4
#     3.01284341153969E-11,  # Scenario 5
#     3.12566420663081E-11,  # Scenario 6
#     2.44722987552953E-11,  # Scenario 7
#     2.00886286074554E-11,  # Scenario 8
#     2.87949683295412E-11,  # Scenario 9
#     2.78124147619842E-11   # Scenario 10
# ])
# adjusting
We = np.array([
    1.107533022806718e-11,  # Scenario 1
    2.998536727832599e-11,   # Scenario 2
    2.441518883756004e-11,   # Scenario 3
    1.126487136125143e-11,  # Scenario 4
    3.436564097908541e-11,  # Scenario 5
    3.11986239282366e-11,  # Scenario 6
    2.014172856588652e-11,  # Scenario 7
    1.534054057620717e-11,  # Scenario 8
    2.924182091720404e-11,  # Scenario 9
    3.222388455288184e-11   # Scenario 10
])
# old shit
# We = np.array([
#     1.045617431394952e-11,  # Scenario 1: C1
#     3.748463916134838e-11,   # Scenario 2: C2
#     2.710237624318629e-11,   # Scenario 3: C3
#     1.23302350017277e-11,   # Scenario 4: C4
#     4.46041974774413E-11,    # Scenario 5: C5
#     3.17336413456305E-11,   # Scenario 6: C6
#     2.00402189519157E-11,   # Scenario 7: C7
#     1.55698207357283E-11,   # Scenario 8: C8
#     2.9353760140396E-11,   # Scenario 9: C9
#     4.23999282033767E-11    # Scenario 10: C10
# ])

# Ensure We is a column vector
We = We.reshape((10, 1))

# Check if M is invertible by computing its determinant
det_M = np.linalg.det(M)
print(f"Determinant of M: {det_M}")

if det_M == 0:
    print("Matrix M is singular and cannot be inverted. Consider revising your simulation scenarios.")
else:
    # Solve for C: C = M^{-1} * (2 * We)
    C = np.linalg.inv(M).dot(2 * We)

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
