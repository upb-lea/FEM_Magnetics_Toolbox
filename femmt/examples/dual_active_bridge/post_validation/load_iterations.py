from Functions import store_as_npy_in_directory
import numpy as np
from matplotlib import pyplot as plt
from operator import itemgetter

plt.figure(figsize=(7, 2.5))
loc = "local_corrected_experimental"
# loc = "local_corrected_extended"
# loc = "local_corrected"
# loc = "local"

for frequency in [200000, 225000, 250000, 275000, 300000]:

    pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/{loc}/results_{frequency}.npy"

    results = np.load(pathD, allow_pickle=True)
    print(len(results))
    print(results)



    # Loss Plots
    total_losses_sorted = [item['total_losses'] for item in results]
    total_losses_sorted = total_losses_sorted[0:10]



    plt.scatter(np.arange(1, len(total_losses_sorted)+1), total_losses_sorted, marker="x", label=f"global grid at {int(frequency/1000)} kHz")




for run in range(1, 2):
    for frequency in [200000]:


        pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/{loc}_{run}/Result_FEM_parameters_{frequency}.npy"
        path_reluctance_parameters = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/{loc}_{run}/reluctance_parameters_{frequency}.npy"
        path_valid_reluctance_parameters = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/{loc}_{run}/valid_reluctance_parameters_{frequency}.npy"
        path_reluctance_parameters_hyst_good = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/{loc}_{run}/reluctance_parameters_hyst_good{frequency}.npy"


        results = np.load(pathD, allow_pickle=True)
        reluctance_parameters = np.load(path_reluctance_parameters, allow_pickle=True)
        valid_reluctance_parameters = np.load(path_valid_reluctance_parameters, allow_pickle=True)
        reluctance_parameters_hyst_good = np.load(path_reluctance_parameters_hyst_good, allow_pickle=True)

        print(f"{len(reluctance_parameters)=}")
        print(f"{valid_reluctance_parameters[0:200]=}")
        # for par in reluctance_parameters:
        #     # if par['p_hyst_nom'] < 8:
        #     #     print(par)
        #
        #     if (par['midpoint'] == 30) and (par["N"] == np.array([[29, 7], [0, 6]])).all():
        #         print(par)

        print(f"{len(valid_reluctance_parameters)=}")
        print(f"{len(reluctance_parameters_hyst_good)=}")
        print(f"{len(results)=}")
        print(results)

        # print(results)

        total_losses_sorted = [item['total_losses'] for item in results]

        # print(total_losses_sorted)

        plot_total_losses_sorted = total_losses_sorted[0:10]
        # plot_total_losses_sorted = total_losses_sorted

        # plt.scatter(np.arange(1, len(plot_total_losses_sorted)+1), plot_total_losses_sorted, color="tab:blue", marker=".", label=f"local grid {run} at {int(frequency/1000)} kHz")

    plt.ylabel(r"$P_\mathrm{total}$ / W")
    plt.xlabel("Sorted Parameter Vectors")
    # plt.legend(loc='upper center', bbox_to_anchor=(1.1, 1.0),
    #           ncol=1, fancybox=True, shadow=True)

    plt.legend(loc='center', bbox_to_anchor=(1.2, 0.5))
plt.grid()
# plt.ylim(15.5, 26.5)

# plt.annotate("^Lab Prototype", (1.8, 16))

plt.savefig(f"C:/Users/tillp/sciebo/Exchange Till/04_Documentation/iterations.png", dpi=1000,
            bbox_inches="tight")
plt.show()
