from Functions import store_as_npy_in_directory
import numpy as np
from matplotlib import pyplot as plt
from operator import itemgetter

# pathA = 'C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/Reluctance_Model/parameter_results.npy'
# pathB = 'C:/Users/tillp/sciebo/Exchange_DAB/current_shapes.npy'
# pathC = 'C:/Users/tillp/Downloads/current_shapes.npy'
result_directory = "C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local"


hyst_1st_FEM = []
hyst_1st_ana = []
hyst_ana = []

# for frequency in [200000, 225000, 250000, 275000, 300000]:
for frequency in [200000]:

    pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/Result_FEM_parameters_{frequency}.npy"
    pathE = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/valid_reluctance_parameters_{frequency}.npy"


    results = np.load(pathD, allow_pickle=True)
    valid_reluctance_parameters = np.load(pathE, allow_pickle=True)

    print(f"{len(results)=}")
    print(f"{len(valid_reluctance_parameters)=}")
    results = [item for item in results if ("FEM_results" not in item)]

    # # Add total losses to each dictionary
    # results = [dict(item, **{'total_losses': item["p_hyst_nom"] +
    #                                          np.sqrt(np.mean(np.array(item["j2H"])**2)) +
    #                                          np.sqrt(np.mean(np.array(item["j2F"])**2))}) for item in results]
    # Add total losses to each dictionary
    results = [dict(item, **{'total_losses': item["p_hyst_nom"] +
                                             np.sum(item["j2H"]) +
                                             np.sum(item["j2F"])}) for item in results]
    # Sort with ascending total losses
    results = sorted(results, key=itemgetter('total_losses'))

    # store_as_npy_in_directory(result_directory, f"results_{frequency}", results)

    # get nominal Hysteresis Losses
    for single_result in results:
        hyst_1st_FEM.append(single_result["p_hyst"][0])
        hyst_1st_ana.append(single_result["p_hyst_nom_1st"])
        hyst_ana.append(single_result["p_hyst_nom"])


    print(f"{len(results)=}")
    print(f"{results=}")

    # Loss Plots
    total_losses_sorted = [item['total_losses'] for item in results]
    plt.scatter(np.arange(1, len(total_losses_sorted)+1), total_losses_sorted, label=f"{int(frequency/1000)} kHz")


plt.ylabel("Total Losses in W")
plt.xlabel("Sorted Parameter Runs")
plt.legend()
plt.grid()
plt.show()

# Analysis

mean_hyst_1st_FEM = np.mean(hyst_1st_FEM)
mean_hyst_1st_ana = np.mean(hyst_1st_ana)

hyst_1st_FEM = np.array(hyst_1st_FEM)
hyst_1st_ana = np.array(hyst_1st_ana)

dist_mean = np.mean((hyst_1st_FEM-hyst_1st_ana)/hyst_1st_FEM)



print(mean_hyst_1st_ana, mean_hyst_1st_FEM, dist_mean)

plt.figure(figsize=(6, 3.5))
plt.plot(hyst_1st_FEM, label="harmonic excitation, FEM")
plt.plot(hyst_1st_ana, label="harmonic excitation, analytical")
plt.plot(hyst_ana, label=r"time-domain fluxes, analytical")
plt.ylabel(r"$P_{\rm hyst}  / \rm W$")
plt.xlabel("Sorted Parameter Runs")
plt.legend()
plt.grid()
plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/compare_hyst.pdf",
            bbox_inches="tight")
plt.show()