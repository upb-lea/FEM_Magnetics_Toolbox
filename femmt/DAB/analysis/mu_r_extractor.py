import femmt
from femmt import MagneticComponent
import numpy as np
import itertools
import re
import csv
from matplotlib import pyplot as plt

plt.figure(figsize=(5, 2.5))

directory = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/data/Loss_Data_PHD_Keuck/mu_r_Plot"
temperatures = [60, 80, 100]
frequency_in_kHz = 300

for temperature in temperatures:
    file_name = f"mu_r_{frequency_in_kHz}kHz_N95_{temperature}C.txt"

    # Get raw data from measurements
    with open(directory + "/" + file_name, newline='\n') as file:
        raw_data = [re.split(r'\t+', line[0]) for line in list(csv.reader(file))]

    raw_data = np.array([[float(i) for i in j] for j in raw_data])

    data = []
    for dat in raw_data:
        if (0.03 < dat[0] < 0.25):
            data.append(dat)

    data = np.array(data)
    print(data[:40, :])

    plt.plot(1000*data[:, 0], data[:, 1], "--", label=f"{temperature} °C")

plt.axhline(y=3000, color='k', linestyle='-', label=r"inital $\mu_\mathrm{r}$")
plt.ylabel(r"$\mu_\mathrm{r}  /  \mu_0$")
plt.xlabel("$B$ / mT")
# plt.ylim(2600, 3600)
plt.legend(ncol=2)
plt.grid()
plt.savefig("C:/Users/tillp/sciebo/Exchange Till/04_Documentation/Core loss/mu_r.pdf", bbox_inches="tight")
plt.show()
