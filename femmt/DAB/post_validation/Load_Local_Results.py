from femmt_functions import store_as_npy_in_directory
import numpy as np
from matplotlib import pyplot as plt
from operator import itemgetter

# pathA = 'C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/Reluctance_Model/parameter_results.npy'
# pathB = 'C:/Users/tillp/sciebo/Exchange_DAB/current_shapes.npy'
# pathC = 'C:/Users/tillp/Downloads/current_shapes.npy'
result_directory = "C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local"


for frequency in [200000]:

    pathD = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/local/Result_FEM_parameters_{frequency}.npy"

    results = np.load(pathD, allow_pickle=True)

    print(f"{len(results)=}")
    results = [item for item in results if ("FEM_results" not in item)]

    # Add total losses to each dictionary
    results = [dict(item, **{'total_losses': item["p_hyst_nom"] +
                                             np.sum(item["j2H"]) +
                                             np.sum(item["j2F"])}) for item in results]
    # Sort with ascending total losses
    results = sorted(results, key=itemgetter('total_losses'))

    # store_as_npy_in_directory(result_directory, f"results_{frequency}", results)


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
