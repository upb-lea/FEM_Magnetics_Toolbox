import numpy as np
# from femmt_functions import NbrStrands
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import pathlib


full_path= 'log/litz_eva/fillfactors8'

layers = [9]  # [4, 5, 6, 9, 10, 11, 12]
radix = [0.05e-3]  # [0.01e-3, 0.05e-3, 0.1e-3]
fillfactors = [0.4, 0.5, 0.6, 0.7]


df_femm = pd.read_json(full_path+'Pv_FEMM_fill_factor.json', convert_axes=False)
df_femmt = pd.read_json(full_path+'Pv_FEMMT_fill_factor.json', convert_axes=False)

df_Pv_femm = df_femm.set_index('Frequencies')
df_Pv_femmt = df_femmt.set_index('Frequencies')

# df_Pv_femm = pd.read_json('o_Pv_FEMM_fill_factor.json', numpy=True, convert_axes=False)
# df_Pv_femmt = pd.read_json('o_Pv_FEMMT_fill_factor.json', numpy=True, convert_axes=False)
# df_Pv_femm = df_Pv_femm[['10, 5e-05, 0.57', '11, 5e-05, 0.69', '4, 0.0001, 0.42', '5, 0.0001, 0.63']]
# df_Pv_femmt = df_Pv_femmt[['10, 5e-05, 0.57', '11, 5e-05, 0.69', '4, 0.0001, 0.42', '5, 0.0001, 0.63']]

# Correct 0Hz error in femm
print(df_Pv_femm)
print(df_Pv_femmt)
print(df_Pv_femm.iloc[:, 0])
df_Pv_femm.iloc[0] = df_Pv_femm.iloc[0] * 2

frequencies = [0, 50, 100, 150, 200, 250]

plt.figure(figsize=(4, 2.5))

fs = 11


for ff, col in [(0.4, 'g'), (0.5, 'b'), (0.6, 'r'), (0.7, 'y')]:

    plt.plot(frequencies, np.array(list(df_Pv_femm[f"femm, 9, 5e-05, {ff}"])), 'o--', color=col, label=f"{int(ff*100)}% (FEMM)")
    plt.plot(frequencies, np.array(list(df_Pv_femmt[f"onelab, 9, 5e-05, {ff}"])), 'x-', color=col, label=f"{int(ff*100)}% (FEMMT)")

# ax = df_Pv_femm.plot(marker="p")
# df_Pv_femmt.plot(ax=ax, marker="o")
# ax.text(25000, 1, 'FEMM', color='r', ha='right', rotation=0, wrap=True)
# ax.text(25000, 0.5, 'FEMMT', color='g', ha='right', rotation=0, wrap=True)

plt.xlabel(r"$f$ / kHz", fontsize=fs)
plt.xticks(fontsize=fs)
plt.ylabel(r"$P_{\mathrm{winding}}$ / W", fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid()
plt.legend(title=r"Fill Factors:", bbox_to_anchor=(0.5, 1.3), loc='center', ncol=2)
# plt.legend(title=r"Fill Factors:")
plt.savefig("Losses_ff.pdf", bbox_inches="tight")
plt.show()



# Error plotting
error = pd.DataFrame([], columns=[])
for col in range(0, len(df_Pv_femmt.columns)):
    error.insert(loc=0, column=f"FF{fillfactors[-col-1]}", value=(df_Pv_femmt.iloc[:, col] - df_Pv_femm.iloc[:, col]).div(df_Pv_femmt.iloc[:, col]))

print(error)
print(list(error["FF0.4"]))


plt.figure(figsize=(4, 4))


for ff, col in [(0.4, 'g'), (0.5, 'b'), (0.6, 'r'), (0.7, 'y')]:

    plt.plot(frequencies, np.array(list(error[f"FF{ff}"])), 'o-', color=col, label=f"{int(ff*100)}%")

# after plotting the data, format the labels
current_values = plt.gca().get_yticks()
plt.gca().set_yticklabels([f'{np.round(x*100, 3)}%'.format(x) for x in current_values])

plt.xlabel(r"$f$ / kHz", fontsize=fs)
plt.xticks(fontsize=fs)
# plt.ylabel(r"$\frac{P_{cond, FEMMT} - P_{cond, FEMM}}{P_{cond, FEMMT}}$", fontsize=14)
plt.ylabel(r"Relative Deviation", fontsize=fs)
plt.yticks(fontsize=fs)

plt.grid()
# plt.legend(title=r"Fill Factors:", title_fontsize='medium')
plt.legend(title=r"Fill Factors:")
plt.savefig("error_ff.pdf", bbox_inches="tight")
plt.show()


# error.plot(title=r"$\frac{P_{v,onelab}-P_{v,femm}}{P_{v,onelab}}$", )

# print("Error: \n", error)
# print("FEMM: \n", df_Pv_femm)
# print("FEMMT: \n", df_Pv_femmt)

# Single Plots
#df_Pv_femm.plot(title="FEMM")
#df_Pv_femmt.plot(title="FEMMT")


# Two Plots in one
