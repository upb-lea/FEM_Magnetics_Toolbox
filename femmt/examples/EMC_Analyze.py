# import pandas as pd
# import numpy as np
#
#
def read_df_from_comsol_3d_plot(link2file, quantity_name: str):
    df = pd.read_csv(link2file, sep='\s+', comment="%", skiprows=8, header=None)
    df.columns = ["x", "y", "z", quantity_name]
    return df

# link2file = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\EMC_Othman\EMC_Othman\Oa1.DAT"
# link2file = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\EMC_Othman\EMC_Othman\Oa2.DAT"

import pandas as pd
import matplotlib.pyplot as plt

def read_df_from_spectrum_dat_file(filepath, quantity_name: str = "Signal Level (dBμV)"):
    with open(filepath, 'r', encoding='latin1') as f:  # <--- changed here
        lines = f.readlines()

    # Find start of the Values section
    start_index = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Values"):
            start_index = i + 1
            break

    if start_index is None:
        raise ValueError("Could not find 'Values' section in file.")

    # Parse data
    data_lines = lines[start_index:]
    data = [list(map(float, line.strip().split(";")[:2])) for line in data_lines if ";" in line]
    df = pd.DataFrame(data, columns=["Frequency (Hz)", quantity_name])
    return df

# Load data from both files
link1 = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\EMC_Othman\EMC_Othman\Oa1.DAT"
link2 = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\EMC_Othman\EMC_Othman\Oa2.DAT"

df1 = read_df_from_spectrum_dat_file(link1)
df2 = read_df_from_spectrum_dat_file(link2)

# Plot
plt.figure(figsize=(12, 6))
plt.plot(df1["Frequency (Hz)"], df1["Signal Level (dBμV)"], label='Oa1.DAT', linewidth=1)
plt.plot(df2["Frequency (Hz)"], df2["Signal Level (dBμV)"], label='Oa2.DAT', linewidth=1)
plt.title("Signal Level vs Frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Signal Level (dBμV)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# df1 = read_df_from_spectrum_dat_file(link1)
# df2 = read_df_from_spectrum_dat_file(link2)
#
# # Convert Hz to MHz for better readability
# df1["Frequency (MHz)"] = df1["Frequency (Hz)"] / 1e6
# df2["Frequency (MHz)"] = df2["Frequency (Hz)"] / 1e6
#
# # Plotting with adjusted style
# plt.style.use("dark_background")
# plt.figure(figsize=(14, 7))
#
# plt.plot(df1["Frequency (MHz)"], df1["Signal Level (dBμV)"], color='yellow', linewidth=1.2, label='Oa1.DAT')
# plt.plot(df2["Frequency (MHz)"], df2["Signal Level (dBμV)"], color='orange', linewidth=1.2, label='Oa2.DAT')
#
# # Red horizontal line at 73 dBμV
# plt.axhline(y=73, color='red', linestyle='--', linewidth=1)
# plt.text(df1["Frequency (MHz)"].min(), 73 + 1, '73 dBμV', color='red', fontsize=10)
#
# # Grid + styling
# plt.title("Signal Level vs Frequency", fontsize=14)
# plt.xlabel("Frequency (MHz)")
# plt.ylabel("Signal Level (dBμV)")
# plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
# plt.legend()
# plt.tight_layout()
# plt.show()

def read_df_from_spectrum_dat_file_2(filepath, quantity_name: str = "Signal Level (dBμV)"):
    with open(filepath, "r", encoding="latin1") as f:
        lines = f.readlines()

    # locate “Values” section
    for i, line in enumerate(lines):
        if line.strip().startswith("Values"):
            start_index = i + 1
            break
    else:
        raise ValueError("Could not find 'Values' section in file.")

    # read frequency / level pairs
    data = [
        list(map(float, line.strip().split(";")[:2]))
        for line in lines[start_index:]
        if ";" in line
    ]
    return pd.DataFrame(data, columns=["Frequency (Hz)", quantity_name])

# load data
link1 = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\EMC_Othman\EMC_Othman\Oa1.DAT"
link2 = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\EMC_Othman\EMC_Othman\Oa2.DAT"

df1 = read_df_from_spectrum_dat_file_2(link1)
df2 = read_df_from_spectrum_dat_file_2(link2)

# plot with log-scale frequency axis
plt.figure(figsize=(12, 6))
plt.semilogx(df1["Frequency (Hz)"], df1["Signal Level (dBμV)"], label="Transformer A", linewidth=1)
plt.semilogx(df2["Frequency (Hz)"], df2["Signal Level (dBμV)"], label="Transformer B", linewidth=1)
plt.title("Signal Level vs Frequency")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Signal Level (dBµV)")
plt.grid(True, which="both", ls="--")
plt.legend()
plt.tight_layout()
out_path = r"C:\Users\uthmn\OneDrive - Universität Paderborn\Drive D\Paderborn Uni\Paderborn Uni\Master Thesis\Pictures\masterthesis\Till_transformer\emc_spectrum_transformers.pdf"
plt.savefig(out_path, format="pdf")
plt.show()


