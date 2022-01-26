import numpy as np

pathE = f"C:/Users/tillp/OneDrive/Documents/GitHub/FEM_Magnetics_Toolbox/femmt/MA/final/reluctance_parameters_200000.npy"
init = np.load(pathE, allow_pickle=True)
print(len(init))
print(init)
