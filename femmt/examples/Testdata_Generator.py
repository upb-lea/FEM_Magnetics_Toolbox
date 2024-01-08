import femmt as fmt
import shutil
from DESSCA import dessca_model
import numpy as np
import traceback

class SharedWindingNumber:
    Value = 100

    def __init__(self):
        pass

    def clear(self):
        self.Value = 1

    def increase_shared_value(self):
        self.Value += 1
    def get_shared_value(self):
        return self.Value


class Model_Values:

    Werte = None

    my_dessca_instance = dessca_model(box_constraints=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1],
                                                       [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]],
                                       state_names=["Frequenz","Strom", "Temperatur", "Kerndurchmesser",
                                                    "FensterHöhe", "Fenster-Breite", "Luftspalthöhe",
                                                    "Luftspaltanzahl", "Kernisolation-Links", "Leiterradius",
                                                    "Wickelfensterspaltung", "Leiteranordnung"],
                                       bandwidth=0.1,
                                       render_online= False,
                                       )

    next_sample_suggest = my_dessca_instance.update_and_sample()

    def __init__(self):
        pass

    def next(self):
        self.next_sample_suggest = self.my_dessca_instance.update_and_sample(np.transpose([self.next_sample_suggest]))
        self.Werte = self.next_sample_suggest


class ChangeFolder:
    simulation_number = 0

    def __init__(self):
        pass

    def copy_result_to_new_folder(self):
        # Pfad zur Result Datei
        origin_path = "C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\parallel\\default\\model_0\\results\\log_electro_magnetic.json"

        # Neuer Pfad und Dateiname
        new_path = f"C:\\Users\\samet\\test\\nach\\results\\Results_{self.simulation_number}.json"

        shutil.copy(origin_path, new_path)

        origin_path = "C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\parallel\\default\\model_0\\results\\log_coordinates_description.json"

        new_path = f"C:\\Users\\samet\\test\\nach\\coordinates\\Coordinates_{self.simulation_number}.json"

        shutil.copy(origin_path, new_path)

        self.simulation_number += 1


def denormalize(normalized_value, min, max):
    denormalized_value = round((normalized_value * (max-min))+min)
    return denormalized_value


def denormalize_milli(normalized_value, min, max):
    denormalized_value = round((normalized_value * (max-min))+min,5)
    return denormalized_value


def frequency():
    frequency_value = denormalize(Model_Values.next_sample_suggest[0], 100000, 500000)

    return frequency_value


def current():
    current_value = [denormalize(Model_Values.next_sample_suggest[1], 1, 50)]
    return current_value


def temperature():
    temperature_value = denormalize(Model_Values.next_sample_suggest[2], 30, 100)

    return temperature_value


def corediameter():
    corediameter_value = denormalize_milli(Model_Values.next_sample_suggest[3], 7e-3, 26e-3)

    return corediameter_value


def coreheight():
    corediameter_value = denormalize_milli(Model_Values.next_sample_suggest[3], 7e-3, 26e-3)
    windowheight_value = denormalize_milli(Model_Values.next_sample_suggest[4], 9.5e-3, 35.9e-3)
    coreheight_value = windowheight_value + corediameter_value / 2

    return coreheight_value


def windowheight():
    windowheight_value = denormalize_milli(Model_Values.next_sample_suggest[4], 9.5e-3, 35.9e-3)
    print (windowheight_value)
    return windowheight_value


def windowwidth():
    windowwidth_value = denormalize_milli(Model_Values.next_sample_suggest[5], 3.7e-3, 19.5e-3)

    return windowwidth_value


def airgapnumber():
    airgapnumber = denormalize_milli(Model_Values.next_sample_suggest[7], 0, 3)

    if 0 <= airgapnumber < 1:
        airgapnumber = 1

    if 1 <= airgapnumber < 2:
        airgapnumber = 2

    if 2 <= airgapnumber <= 3:
        airgapnumber = 3
    airgapnumber = 1
    return airgapnumber

def airgapposition(airgap):
    airgapnumber = denormalize_milli(Model_Values.next_sample_suggest[7], 0, 3)
    airgapposition_value = None

    if 0 <= airgapnumber < 1:
        airgapposition_value = 50

    if 1 <= airgapnumber < 2:
        if airgap == 1:
            airgapposition_value = 33.33

        if airgap == 2:
            airgapposition_value = 66.66

    if 2 <= airgapnumber <= 3:
        if airgap == 1:
            airgapposition_value = 25

        if airgap == 2:
            airgapposition_value = 50

        if airgap == 3:
            airgapposition_value = 75

    return airgapposition_value


def airgapheight():
    Samplevalue = Model_Values.next_sample_suggest[6]

    if Samplevalue <= 0.5:
        airgapheight_value = 0.0005

    if Samplevalue > 0.5:
        airgapheight_value = 0.00125

    return airgapheight_value


def innerwindinginsulation():
    innerwindinginsulation_value = 0.00095

    return innerwindinginsulation_value


def coreinsulation():
    coreinsulation_value = [0.001, 0.001, denormalize_milli(Model_Values.next_sample_suggest[8], 0.0009, 0.001), 0.001]

    return coreinsulation_value


def conductorradius():
    conductorradius_value = denormalize_milli(Model_Values.next_sample_suggest[9], 0.00055, 0.0023)

    return conductorradius_value


# def windingwindowsplit():
#     split = get_value("WindingWindowSplit")
#     windingwindowsplit_value = None
#
#     if split == 1:
#         windingwindowsplit_value = fmt.WindingWindowSplit.NoSplit
#
#     if split == 2:
#         windingwindowsplit_value = fmt.WindingWindowSplit.HorizontalAndVerticalSplit
#
#     return windingwindowsplit_value


def conductorarrangement():
    arrangement = denormalize_milli(Model_Values.next_sample_suggest[11], 0, 3)
    conductorarrangement_value = None

    if 0 <= arrangement < 1:
        conductorarrangement_value = fmt.ConductorArrangement.Square

    if 1 <= arrangement < 2:
        conductorarrangement_value = fmt.ConductorArrangement.SquareFullWidth

    if 2 <= arrangement < 3:
        conductorarrangement_value = fmt.ConductorArrangement.Hexagonal

    return conductorarrangement_value


if __name__ == "__main__":
    Combination = 0
    while Combination < 1000000:
        print("Combination:", Combination)
        # exec(open("C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\hpc_examples.py").read())
        try:
            exec(open("C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\hpc_examples.py").read())
        except Exception as e:
            traceback.print_exc()

        Combination = Combination + 1
