import femmt
from femmt import MagneticComponent
import numpy as np
import itertools
from matplotlib import pyplot as plt

pathA = 'C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/Reluctance_Model/parameter_results.npy'
pathB = 'C:/Users/tillp/sciebo/Exchange_DAB/current_shapes.npy'
pathC = 'C:/Users/tillp/Downloads/current_shapes.npy'

for frequency in [200000, 225000, 250000, 275000, 300000]:

    pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/Result_FEM_parameters_{frequency}.npy"

    A = np.load(pathD, allow_pickle=True)

    # print(A[0]['window_w'])
    print(len(A))
    # print(A)

    hyst = [res["p_hyst_nom"] for res in A]

    print(hyst)
    print(len(hyst))

    hyst_sorted = sorted(hyst)
    print(f"{hyst_sorted=}")

    # max_ok = int(len(hyst_sorted)/10)
    # hyst_good = hyst_sorted[0:max_ok]
    epsilon = 1.25
    hyst_good = [item for item in hyst_sorted if item < hyst_sorted[0]*epsilon]  #TODO: do before FEM simulation, if sth like 150% of minimum ...



    print(hyst_good)
    print(len(hyst_good))

    plt.plot(hyst_good, label=frequency)

plt.legend()
plt.show()