import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pathlib
import fileinput
import csv

# path = str(pathlib.Path(__file__).parent.absolute())
path = "C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/Strands_Coefficents/coeff/"

layers = [4]
fillfactors = [0.4, 0.45, 0.5, 0.55, 0.6, 0.63, 0.71, 0.76]
X_choose = 1.0

results = [[], [], [], []]

for FF in fillfactors:
    for n_layers in layers:
        # Formatting stuff
        files = [path + f"/pB_RS_la{FF}_{n_layers}layer.dat",
                 path + f"/pI_RS_la{FF}_{n_layers}layer.dat",
                 path + f"/qB_RS_la{FF}_{n_layers}layer.dat",
                 path + f"/qI_RS_la{FF}_{n_layers}layer.dat"]

        for i in range(0, 4):
            with open(files[i], newline='\n') as f:
                count = 0
                for line in f:
                    count += 1
                    line = line.strip()
                    words = line.split(sep=' ')
                    words = [float(i) for i in words]
                    results[i].append(words)
                    print(count)

results = np.array(results)
vals = [[], [], [], []]
for num in range(0, 4):
    for X, val in results[num]:
        if X == X_choose:
            vals[num].append(val)

# As numpy array
vals = np.asarray(vals)



# Normalized on average


print(vals)
plt.figure(figsize=(2.5, 3))

plt.plot(fillfactors, vals[0] / np.mean(vals[0]), 'o-', label=f'pB')
plt.plot(fillfactors, vals[1] / np.mean(vals[1]), 'o-', label=f'pI')
plt.plot(fillfactors, vals[2] / np.mean(vals[2]), 'o-', label=f'qB')
plt.plot(fillfactors, vals[3] / np.mean(vals[3]), 'o-', label=f'qI')

# after plotting the data, format the labels
current_values = plt.gca().get_yticks()
plt.gca().set_yticklabels([f'{np.round(x*100-100, 3)}%'.format(x) for x in current_values])

from matplotlib.pyplot import figure

# plt.figure(figsize=(6, 8), dpi=80)
# plt.ylim([0.999, 1.001])
plt.xlabel('Fill Factors')
plt.ylabel(r"Normalized on averages")
# plt.title('Relative deviation of strand parameters \n normalized on their averages')
plt.legend()
plt.grid()
plt.savefig("C:/Users/tillp/sciebo/Exchange Till/04_Documentation/Inkscape/Litz wires/strands_ff_all.pdf", bbox_inches="tight")
plt.show()