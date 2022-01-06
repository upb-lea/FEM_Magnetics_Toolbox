from femmt_functions import store_as_npy_in_directory
import numpy as np
from matplotlib import pyplot as plt
from operator import itemgetter

plt.figure(figsize=(6, 4))

for frequency in [200000, 225000, 250000, 275000, 300000]:

    pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local/results_{frequency}.npy"

    results = np.load(pathD, allow_pickle=True)
    print(len(results))



    # Loss Plots
    total_losses_sorted = [item['total_losses'] for item in results]
    plt.scatter(np.arange(1, len(total_losses_sorted)+1), total_losses_sorted, label=f"Gridsearch at {int(frequency/1000)} kHz")


for run in range(1, 4):
    for frequency in [200000]:

        pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_{run}/Result_FEM_parameters_{frequency}.npy"
        path_reluctance_parameters = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_{run}/reluctance_parameters_{frequency}.npy"
        path_valid_reluctance_parameters = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_{run}/valid_reluctance_parameters_{frequency}.npy"
        path_reluctance_parameters_hyst_good = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local_{run}/reluctance_parameters_hyst_good{frequency}.npy"


        results = np.load(pathD, allow_pickle=True)
        reluctance_parameters = np.load(path_reluctance_parameters, allow_pickle=True)
        valid_reluctance_parameters = np.load(path_valid_reluctance_parameters, allow_pickle=True)
        reluctance_parameters_hyst_good = np.load(path_reluctance_parameters_hyst_good, allow_pickle=True)

        print(f"{len(reluctance_parameters)=}")
        print(f"{len(valid_reluctance_parameters)=}")
        print(f"{len(reluctance_parameters_hyst_good)=}")
        print(f"{len(results)=}")

        # print(results)

        total_losses_sorted = [item['total_losses'] for item in results]

        # print(total_losses_sorted)

        plt.scatter(np.arange(1, len(total_losses_sorted)+1), total_losses_sorted, label=f"local run {run} at {int(frequency/1000)} kHz")

    plt.ylabel(r"$P_\mathrm{total}$ / W")
    plt.xlabel("Sorted Parameter Vectors")
    plt.legend(loc='upper center', bbox_to_anchor=(1.1, 1.0),
              ncol=1, fancybox=True, shadow=True)
    # plt.legend(ncol=2)
plt.grid()
plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/iterations.pdf",
            bbox_inches="tight")
plt.show()
