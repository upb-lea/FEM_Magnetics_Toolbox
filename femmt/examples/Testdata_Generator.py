import femmt as fmt
import shutil
from itertools import product


"""------------------------COMPONENT---------------------------"""
"Component: Transformator = 1 or Inductor = 2"
component = 2

# Frequency_Arr = [100000, 200000]
# Current_Arr = [1, 8, 15]
# Temperature_Arr = [30, 60,80, 100]
# CoreDiameter_Arr = [0.007, 8.8e-3]
# WindowHeight_Arr = [9.5e-3, 16.9e-3]
# WindowWidth_Arr = [3.7e-3, 6.4e-3]
# AirGapHeight_Arr = [0.0005]
# AirGapNumber_Arr = [1, 2, 3]  # Number of AirGaps
# CoreInsulation_LeftCore_Arr = [0.001]  # [0.001, 0.001, CoreInsulation_LeftCore, 0.001]
# ConductorRadius_Arr = [0.00055]
# WindingWindowSplit_Arr = [1]  # NoSplit = 1, HorizontalAndVerticalSplit = 2
# ConductorArrangement_Arr = [1]  # Square = 1, SquareFullWidth = 2, Hexagonal = 3

# Werte der Parameter
Frequency_Arr = [100000, 200000, 300000, 400000, 500000]
Current_Arr = [1, 12.5, 25, 37.5, 50]
Temperature_Arr = [30, 60,80, 100]
CoreDiameter_Arr = [7e-3, 14.9e-3, 20e-3, 26e-3]
WindowHeight_Arr = [9.5e-3, 22.9e-3, 35.9e-3]
WindowWidth_Arr = [3.7e-3, 9.7e-3, 19.5e-3]
AirGapHeight_Arr = [0.0005, 0.00125]
AirGapNumber_Arr = [1, 2, 3]  # Number of AirGaps
CoreInsulation_LeftCore_Arr = [0.0009, 0.00095, 0.001]  # [0.001, 0.001, CoreInsulation_LeftCore, 0.001]
ConductorRadius_Arr = [0.00055, 0.0015, 0.0023]
WindingWindowSplit_Arr = [1, 2]  # NoSplit = 1, HorizontalAndVerticalSplit = 2
ConductorArrangement_Arr = [1, 2, 3]  # Square = 1, SquareFullWidth = 2, Hexagonal = 3


# Erstellen aller möglichen Kombinationen
combinations = list(product(Frequency_Arr, Current_Arr, Temperature_Arr,CoreDiameter_Arr, WindowHeight_Arr, WindowWidth_Arr, AirGapHeight_Arr, AirGapNumber_Arr, CoreInsulation_LeftCore_Arr, ConductorRadius_Arr, WindingWindowSplit_Arr, ConductorArrangement_Arr))
print(len(combinations))

# Funktion gibt die Position des entsprechenden Parameters in der Kombination zurück
def get_number_for_list(input_string):
    array_mapping = {
        "Frequency": 0,
        "Current": 1,
        "Temperature": 2,
        "CoreDiameter": 3,
        "WindowHeight" : 4,
        "WindowWidth" : 5,
        "AirGapHeight" : 6,
        "AirGapNumber" : 7,
        "CoreInsulation_LeftCore" : 8,
        "ConductorRadius" : 9,
        "WindingWindowSplit" : 10,
        "ConductorArrangement": 11
    }
    default = []

    return array_mapping.get(input_string, default)


# Funktion verschiebt die Ergebnisse in einen anderen Ordner
def copy_result_to_new_folder(simulation_number):
    # Pfad zur Result Datei
    origin_path = "C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\parallel\\default\\model_0\\results\\log_electro_magnetic.json"

    # Neuer Pfad und Dateiname
    new_path = f"C:\\Users\\samet\\test\\nach\\Results_{simulation_number}.json"

    shutil.copy(origin_path, new_path)

# Funktion gibt den Wert des Parameters in der entsprechenden Simulation zurück
def get_value(variable):
    # Parameterkombination für die jeweilige Simulation
    combination = combinations[SharedSimulationNumber.Value]

    # Parameter der gesuchten Variable für die jeweilige Simulation
    value = combination[get_number_for_list(variable)]
    return value


def frequency():
    frequency_value = get_value("Frequency")

    return frequency_value


def current():
    current_value = None

    if component == 1:
        current_value = [get_value("Current"), 10]

    if component == 2:
        current_value = [get_value("Current")]

    return current_value


def temperature():
    temperature_value = get_value("Temperature")

    return temperature_value


def corediameter():
    corediameter_value = get_value("CoreDiameter")
    return corediameter_value


def coreheight():
    corediameter_value = get_value("CoreDiameter")
    windowheight_value = get_value("WindowHeight")
    coreheight_value = windowheight_value + corediameter_value / 2

    return coreheight_value


def windowheight():
    windowheight_value = get_value("WindowHeight")

    return windowheight_value


def windowwidth():
    windowwidth_value = get_value("WindowWidth")

    return windowwidth_value


def airgapposition(airgap):
    airgapnumber = get_value("AirGapNumber")
    airgapposition_value = None

    if airgapnumber == 1:
        airgapposition_value = 50

    if airgapnumber == 2:
        if airgap == 1:
            airgapposition_value = 33.33

        if airgap == 2:
            airgapposition_value = 66.66

    if airgapnumber == 3:
        if airgap == 1:
            airgapposition_value = 25

        if airgap == 2:
            airgapposition_value = 50

        if airgap == 3:
            airgapposition_value = 75

    return airgapposition_value

def airgapheight():
    airgapheight_value = get_value("AirGapHeight")

    return airgapheight_value


def innerwindinginsulation():
    innerwindinginsulation_value = None

    if component == 1:
        innerwindinginsulation_value = [[0.0002, 0.001], [0.001, 0.0002]]

    if component == 2:
        innerwindinginsulation_value = 0.00095

    return innerwindinginsulation_value


def coreinsulation():
    coreinsulation_value = [0.001, 0.001, get_value("CoreInsulation_LeftCore"), 0.001]

    return coreinsulation_value


def conductorradius():
    conductorradius_value = get_value("ConductorRadius")

    return conductorradius_value


def windingwindowsplit():
    split = get_value("WindingWindowSplit")
    windingwindowsplit_value = None

    if split == 1:
        windingwindowsplit_value = fmt.WindingWindowSplit.NoSplit

    if split == 2:
        windingwindowsplit_value = fmt.WindingWindowSplit.HorizontalAndVerticalSplit

    return windingwindowsplit_value


def windingturns(winding):
    windingturns_value = None
    wind_arr = [1,1,2,3]
    if winding == 1:
        windingturns_value = 3

    if winding == 2:
        windingturns_value = 21


    return windingturns_value


def conductorarrangement():
    arrangement = get_value("ConductorArrangement")
    conductorarrangement_value = None

    if arrangement == 1:
        conductorarrangement_value = fmt.ConductorArrangement.Square

    if arrangement == 2:
        conductorarrangement_value = fmt.ConductorArrangement.SquareFullWidth

    if arrangement == 3:
        conductorarrangement_value = fmt.ConductorArrangement.Hexagonal

    return conductorarrangement_value


# Klasse überträgt Programmübergreifend die aktuelle Simulationsnummer
class SharedSimulationNumber:
    Value = -1  # Die gemeinsame Variable wird in der Klasse initialisiert

    def __init__(self):
        # Initialisiere die gemeinsame Variable nur, wenn sie noch nicht initialisiert wurde
        pass

    def set_shared_value(self, value):
        # Diese Methode setzt den Wert der gemeinsamen Variable
        SharedSimulationNumber.Value = value

    def get_shared_value(self):
        # Diese Methode gibt den aktuellen Wert der gemeinsamen Variable zurück
        return SharedSimulationNumber.Value


if __name__ == "__main__":
    Simulation = 0
    while Simulation < len(combinations):
        print(Simulation)
        try:
            exec(open("C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\hpc_examples.py").read())
            copy_result_to_new_folder(Simulation)

        except:
            print("Sorry")

        Simulation = Simulation + 1