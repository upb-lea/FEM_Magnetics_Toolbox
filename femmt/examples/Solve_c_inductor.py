
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
We_case_1 = np.array([
    5.23950306321443e-12,  # Scenario 1
    1.420017752187458e-11, # Scenario 2
    5.05540791441735e-12,# Scenario 3
])
# inductor 84
# We_case_1 = np.array([
#     2.516035519332448e-11,  # Scenario 1
#     1.480500554818899e-11, # Scenario 2
#     1.99302692819917e-11,# Scenario 3
# ])
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
C_AB_42 = -3.91159e-12 + (1.43906e-11 * 1.40098e-11 / (1.43906e-11 + 1.40098e-11))
print("C_42:", C_AB_42)
# C_AB_84 = 3.02856e-11 + (2.00351e-11 * 9.57492e-12 / (2.00351e-11 + 9.57492e-12))
# print("C_84:", C_AB_84)
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
import matplotlib.pyplot as plt
import numpy as np

# Data
parameters = ['C1', 'C2', 'C3', 'C_AB(64)']
vertical_I = [32.35, 22.49, 13.65, 40.85]
vertical_II = [22.45, 23.32, 12.82, 30.72]
vertical_III = [29.74, 23.02, 11.63, 37.47]

# Plot settings
bar_width = 0.2
x = np.arange(len(parameters))

plt.figure(figsize=(10, 5))
plt.bar(x - bar_width, vertical_I, width=bar_width, label='Vertical I')
plt.bar(x, vertical_II, width=bar_width, label='Vertical II')
plt.bar(x + bar_width, vertical_III, width=bar_width, label='Vertical III')

plt.xticks(x, parameters)
plt.ylabel('Capacitance (pF)')
plt.title('Comparison of capacitance values for different vertical arrangements')
plt.legend()
plt.grid(True, axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Save as PDF BEFORE showing the plot
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\Final presentation\vertical_comparison.pdf"
plt.savefig(pdf_path)

plt.show()


# Data
parameters = ['C1', 'C2', 'C3', 'C_AB(64)']
Horizontal_I = [2.28, 31.25, 14.38, 12.13]
Horizontal_II = [0.90, 30.85, 14.08, 10.57]
Horizontal_III = [1.26, 30.96, 12.82, 10.33]

# Plot settings
bar_width = 0.2
x = np.arange(len(parameters))

plt.figure(figsize=(10, 5))
plt.bar(x - bar_width, Horizontal_I, width=bar_width, label='Horizontal I')
plt.bar(x, Horizontal_II, width=bar_width, label='Horizontal II')
plt.bar(x + bar_width, Horizontal_III, width=bar_width, label='Horizontal III')

plt.xticks(x, parameters)
plt.ylabel('Capacitance (pF)')
plt.title('Comparison of capacitance values for different horizontal arrangements')
plt.legend()
plt.grid(True, axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

# Save as PDF BEFORE showing the plot
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\Final presentation\Horizontal_comparison.pdf"
plt.savefig(pdf_path)

plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Capacitance labels
cap_labels = [f'C{i}' for i in range(1, 11)]

# Capacitance values from the image
transformer_A = [-6.65, -4.27, 4.27, 3.97, 9.48, 9.82, 12.14, 11.11, 7.12, 6.03]
transformer_B = [-5.60, -6.07, 6.89, 7.05, 3.61, 3.04, 12.63, 11.57, 17.79, 16.56]

# Plot settings
x = np.arange(len(cap_labels))
bar_width = 0.35

plt.figure(figsize=(12, 5))
plt.bar(x - bar_width/2, transformer_A, width=bar_width, label='Transformer A')
plt.bar(x + bar_width/2, transformer_B, width=bar_width, label='Transformer B')

# Labeling
plt.xticks(x, cap_labels)
plt.ylabel('Capacitance (pF)')
plt.title('Comparison of capacitances for transformer A and B')
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()

# Save as PDF BEFORE showing the plot
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\Final presentation\transformer_a_b.pdf"
plt.savefig(pdf_path)

plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Parameters
labels = ['C1 (pF)', 'C2 (pF)', 'C3 (pF)', 'C_AB (pF)']
single_layer = [-3.84, 16.94, 15.24, 4.18]
two_layers = [32.35, 22.49, 13.65, 40.85]
three_layers = [28.83, 26.16, 14.82, 38.07]

# Bar plot positions
x = np.arange(len(labels))
bar_width = 0.25

# Create plot
plt.figure(figsize=(10, 5))
plt.bar(x - bar_width, single_layer, width=bar_width, label='Single-layer')
plt.bar(x, two_layers, width=bar_width, label='Two-layers')
plt.bar(x + bar_width, three_layers, width=bar_width, label='Three-layers')

# Labeling
plt.xticks(x, labels)
plt.ylabel('Capacitance (pF)')
plt.title('Comparison of capacitance for single, two, and three-layer inductors')
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()

plt.show()

# Labels and data
labels = ['C1 (pF)', 'C2 (pF)', 'C3 (pF)', 'C_AB (pF)', 'Inductance (mH)', 'Losses (W)']
single_layer = [-3.84, 16.94, 15.24, 4.18, 0.503, 6.74]
two_layers = [32.35, 22.49, 13.65, 40.85, 2.012, 28.756]
three_layers = [28.83, 26.16, 14.82, 38.07, 4.27, 67.33]
four_layers = [24.71, 28.39, 16.83, 35.28, 8.05, 123.56]
five_layers = [20.61, 31.01, 19.64, 32.26, 12.58, 198.29]

# X positions
x = np.arange(len(labels))
bar_width = 0.15  # smaller bar width for more groups

# Plot
plt.figure(figsize=(14, 6))
plt.bar(x - 2*bar_width, single_layer, width=bar_width, label='Single-layer')
plt.bar(x - bar_width, two_layers, width=bar_width, label='Two-layers')
plt.bar(x, three_layers, width=bar_width, label='Three-layers')
plt.bar(x + bar_width, four_layers, width=bar_width, label='Four-layers')
plt.bar(x + 2*bar_width, five_layers, width=bar_width, label='Five-layers')

# Labels and formatting
plt.xticks(x, labels, rotation=45)
plt.ylabel('Values')
plt.title('Comparison of capacitance, inductance, and losses across layer configurations')
plt.legend()
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()

# Save to PDF
pdf_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\Final presentation\inductor_five_layers.pdf"
plt.savefig(pdf_path)

plt.show()







