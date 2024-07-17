import femmt as fmt
import os

copper = fmt.pcbvariables.Layers(LayerThickness = 35e-6, LayerNumber = 6, WindingPerLayer = 2)

def copper_check(LayerThickness, LayerNumber, WindingPerLayer):
    if LayerThickness not in [35e-6, 70e-6, 105e-6]:
        raise Exception("Copper Thickness no standard Value.")
    if LayerNumber % 2:
        raise Exception("Enter even Number of Layers")
print(copper)

