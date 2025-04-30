
import numpy as np
import matplotlib.pyplot as plt

# Define the M matrix based on your corrected scenarios

M_case_1 = np.array([
    [1, 1, 0],  # Scenario 1
    [0, 1, 1],  # Scenario 2
    [1, 0, 1]  # Scenario 3
])

M_squared1 = M_case_1 ** 2
# inductor 42
# We_case_1 = np.array([
#     5.23950306321443e-12,  # Scenario 1
#     1.420017752187458e-11, # Scenario 2
#     5.05540791441735e-12,# Scenario 3
# ])
# inductor 84
We_case_1 = np.array([
    2.516035519332448e-11,  # Scenario 1
    1.480500554818899e-11, # Scenario 2
    1.99302692819917e-11,# Scenario 3
])
# inductor 81
# We_case_1 = np.array([
#     2.432750994823512e-11,  # Scenario 1
#     1.812344522204704e-11, # Scenario 2
#     1.751673625287478e-11# Scenario 3
# ])
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

connection_sums = {
    'C_AvsBE': lambda C: C[0] + C[1],
    'C_BvsAE': lambda C: C[0] + C[2],
    'C_ABvsE': lambda C: C[1] + C[2],
    'C_AvsB': lambda C: C[0] + 0.25 * C[1] + 0.25 * C[2],
}

connection_measurement = {
    'C_AvsBE': 1.172807229E-11,
    'C_BvsAE': 1.004441933E-11,
    'C_ABvsE': 2.746247081E-11,
    'C_AvsB': 2.804276221E-12
}

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
    # C_AB_42 = -3.91159e-12 + (1.43906e-11 * 1.40098e-11 / (1.43906e-11 + 1.40098e-11))
    # print("C_42:", C_AB_42)
C_AB_84 = 3.02856e-11 + (2.00351e-11 * 9.57492e-12 / (2.00351e-11 + 9.57492e-12))
print("C_84:", C_AB_84)
    # C_AB_81= 2.37208e-11 + (2.49342e-11 * 1.13127e-11 / (2.49342e-11 + 1.13127e-11))
    # print("C_81:", C_AB_81)
    #
    # # c_ab = (7/5) * 2.096884260266433e-12
    # # print(c_ab)
    # # error = ((3.1872335334713594e-12 - 2.096884260266433e-12) / 3.1872335334713594e-12) * 100
    # # print(error)
    # # error_2 = 3.1872335334713594e-12 / 2.096884260266433e-12
    # # print(error_2)
    # c_tt = 1.740398743702928e-12
    # c_ab_84 = 1.83 * 1.740398743702928e-12
    #
    # c_ab_sim = 1.838169938418392e-11 * 2
    # print(c_ab_84)
    # print(c_ab_sim)
    #
    # C_AB_Othman = (-3.91159e-12 + 0.25 * 1.43906e-11 + 0.25 * 1.40098e-11)
    # print("C_Othman:", C_AB_Othman)
    #
    # print("A vs BE:", -3.91159e-12 + 1.43906e-11)
    # print("B vs AE:", -3.91159e-12 + 1.40098e-11)
    # print("AB vs E:", 1.43906e-11 + 1.40098e-11)



# Compute capacitances (if not already done)
C_case1 = np.linalg.inv(M_case_1 ** 2).dot(2 * We_case_1).flatten()

# Compute simulated and measured capacitance sums in pF
calculated_sums_pf = [connection_sums[key](C_case1) * 1e12 for key in connection_sums]
measured_sums_pf = [connection_measurement[key] * 1e12 for key in connection_sums]

# Clean x-axis labels
connection_labels_clean = [
    "A vs BE", "B vs AE", "AB vs E", "A vs B"
]

# Plot
plt.figure(figsize=(10, 6))
x_indices = np.arange(len(connection_labels_clean))

plt.scatter(x_indices, calculated_sums_pf, label="Simulated Capacitance", color='blue', marker='o')
plt.scatter(x_indices, measured_sums_pf, label="Measured Capacitance", color='red', marker='x')

# Annotate error ratios
for i, (calc, meas) in enumerate(zip(calculated_sums_pf, measured_sums_pf)):
    ratio = calc / meas
    plt.text(i, calc, f"{ratio:.2f}×", fontsize=10, ha='center', va='bottom', color='blue')

# Formatting
plt.xticks(x_indices, connection_labels_clean, rotation=0)
plt.ylabel('Capacitance (pF)')
plt.title('Comparison of Simulated and Measured Capacitance – Inductor')
plt.legend(loc='upper right')
plt.tight_layout()
plt.grid(True)

# Save as PDF
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\inductor_42\inductor_42_error.pdf"
plt.savefig(pdf_path)

plt.show()





