import femmt as fmt
import shutil
from DESSCA import dessca_model
import numpy as np
import traceback
import math
import time
from random import randrange

# # This class is used to store and transmit the number of turns with the hpc_examples.py
# class SharedWindingNumber:
#     Value = 1
#
#     def __init__(self):
#         pass
#
#     def clear(self):
#         self.Value = 1
#
#     def increase_shared_value(self):
#         self.Value += 1
#     def get_shared_value(self):
#         return self.Value

# This class generates the normalized Values for the Current and the Geometry of the Inductor
class Model_Values:

    Werte = None

    my_dessca_instance = dessca_model(box_constraints=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1],
                                                       [0, 1], [0, 1], [0, 1]],
                                       state_names=["Strom", "Kerndurchmesser",
                                                    "FensterHöhe", "Fenster-Breite", "Luftspalthöhe",
                                                    "Kernisolation-Links", "Leiterradius", "Leiteranordnung"],
                                       bandwidth=0.1,
                                       render_online= False,
                                       )

    next_sample_suggest = my_dessca_instance.update_and_sample()

    def __init__(self):
        pass

    def next(self):
        self.next_sample_suggest = self.my_dessca_instance.update_and_sample(np.transpose([self.next_sample_suggest]))
        self.Werte = self.next_sample_suggest

# This
class ChangeFolder:
    simulation_number = 0

    def __init__(self):
        pass

    def copy_result_to_new_folder(self,number_of_models):
         for i in range(number_of_models):
            # Pfad zur Result Datei
            #origin_path = f"C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\parallel\\default\\model_{i}\\results\\log_electro_magnetic.json"
            origin_path = f"/home/upb/stekiner/Bachelorarbeit/FEMMT/FEM_Magnetics_Toolbox/femmt/examples/example_results/parallel/default/model_{i}/results/log_electro_magnetic.json"

            # Neuer Pfad und Dateiname
            #new_path = f"C:\\Users\\samet\\Bachelorarbeit\\Probe\\res\\Results_{self.simulation_number}.json"
            new_path = f"/home/upb/stekiner/Bachelorarbeit/Results/added_log_electro_magnetic/Results_{self.simulation_number}.json"

            shutil.copy(origin_path, new_path)

            # origin_path = f"C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\example_results\\parallel\\default\\model_{i}\\results\\log_electro_magnetic.json"
            origin_path = f"/home/upb/stekiner/Bachelorarbeit/FEMMT/FEM_Magnetics_Toolbox/femmt/examples/example_results/parallel/default/model_{i}/results/log_coordinates_description.json"

            # Neuer Pfad und Dateiname
            # new_path = f"C:\\Users\\samet\\Bachelorarbeit\\Probe\\res\\Results_{self.simulation_number}.json"
            new_path = f"/home/upb/stekiner/Bachelorarbeit/Results/added_log_coordinates_description/Coordinates_{self.simulation_number}.json"

            shutil.copy(origin_path, new_path)

            self.simulation_number += 1


# Function to denormalize the generated Current and Geometry Values and round them to the wanted digits
def denormalize_and_round(normalized_value, min, max, digits):
    denormalized_value = round((normalized_value * (max-min))+min, digits)

    return denormalized_value


def frequencies():
    # Set the suggested frequencies
    frequency_values = [100000, 200000, 300000, 400000, 500000]
    return frequency_values


def temperatures():
    # Set the suggested temperatures
    temperature_values = [30, 60, 80, 100]
    return temperature_values


def current():
    # Denormalize and round the suggested current value
    current_value = denormalize_and_round(Model_Values.next_sample_suggest[0], 0, 0.6, 2)
    return current_value


def corediameter():
    # Denormalize and round the suggested core diameter value
    corediameter_value = denormalize_and_round(Model_Values.next_sample_suggest[1], 7e-3, 26e-3, 4)
    return corediameter_value


def coreheight():
    # Calculate the suggested core height value
    corediameter_value = denormalize_and_round(Model_Values.next_sample_suggest[1], 7e-3, 26e-3,4)
    windowheight_value = denormalize_and_round(Model_Values.next_sample_suggest[2], 9.5e-3, 35.9e-3,4)

    coreheight_value = round(windowheight_value + (corediameter_value / 2), 4)

    return coreheight_value


def windowheight():
    # Denormalize and round the suggested window height value
    windowheight_value = denormalize_and_round(Model_Values.next_sample_suggest[3], 9.5e-3, 35.9e-3,4)
    return windowheight_value


def windowwidth():
    # Denormalize and round the suggested window width value
    windowwidth_value = denormalize_and_round(Model_Values.next_sample_suggest[4], 3.7e-3, 19.5e-3, 4)
    return windowwidth_value


def airgapheight():
    # Denormalize and round the suggested air gap height value
    airgapheight_value = denormalize_and_round(Model_Values.next_sample_suggest[5], 0.0001, 0.0015, 4)
    return airgapheight_value


def coreinsulation():
    coreinsulation_left = denormalize_and_round(Model_Values.next_sample_suggest[6], 0.0005, 0.002, 4)
    coreinsulation_value = [0.001, 0.001, coreinsulation_left, 0.001]

    return coreinsulation_value


def conductorradius():
    conductorradius_value = denormalize_and_round(Model_Values.next_sample_suggest[7], 0.00055, 0.001,5)

    return conductorradius_value


def innerwindinginsulation():
    innerwindinginsulation_value = 0.00095
    return innerwindinginsulation_value


def airgapnumber_and_position():
    number_of_airgaps = randrange(1, 4)

    airgap_positions = {
        1: [50],
        2: [33.33, 66.66],
        3: [25, 50, 75],
    }

    airgap_position = airgap_positions.get(number_of_airgaps, [])

    return number_of_airgaps, airgap_position


def conductorarrangement():
    arrangement = randrange(1, 4)

    conductorarrangements = {
        1: fmt.ConductorArrangement.Square,
        2: fmt.ConductorArrangement.SquareFullWidth,
        3: fmt.ConductorArrangement.Hexagonal,
    }

    conductorarrangement = conductorarrangements.get(arrangement, [])

    return conductorarrangement


def is_magnetic_flux_density_below_limit(turns, core_height, current, core_permeability=3000, limit=3):

    magnetic_field_limit = limit  # Grenzwert für die magnetische Flussdichte in Tesla
    mu_0 = 4 * 3.14159265358979323846e-7

    magnetic_field_strength = (turns * current) / core_height
    magnetic_flux_density = core_permeability * mu_0 * magnetic_field_strength
    print("Turns: ", turns, "Coreheight: ", core_height, "Current: ", current, "H: ", magnetic_field_strength, "B: ", magnetic_flux_density)

    # Überprüfen, ob die magnetische Flussdichte unter dem Grenzwert liegt
    if magnetic_flux_density < magnetic_field_limit:
        return True
    else:
        return False


def check_windings_fit(turns, conductor_radius, insulation_top, insulation_bottom, insulation_left, insulation_right, window_height, window_width, innerwindinginsulation=0.00095):
    # Berücksichtigen der Isolierung am Kern
    effective_height = window_height - (insulation_top + insulation_bottom)
    effective_width = window_width - (insulation_left + insulation_right)

    # Wie viel Leiter passen in der Höhe
    conductors_height = effective_height / ((2 * conductor_radius) + innerwindinginsulation)

    conductors_height = math.floor(conductors_height)

    # Wie viel Leiter passen in der Breite
    conductors_width = effective_width / ((2 * conductor_radius) + innerwindinginsulation)
    conductors_width = math.floor(conductors_width)

    # Wie viele Leiter passen insgesamt
    conductors_total = conductors_width * conductors_height

    #Check ob die Windungen in das Wickelfenster passen
    if turns < conductors_total:
        return True
    if turns > conductors_total:
        return False



if __name__ == "__main__":
    Combination = 0
    while Combination < 1000000:
        print("Combination:", Combination)
        # exec(open("C:\\Users\\samet\\FEMMT\\GIT\\FEM_Magnetics_Toolbox\\femmt\\examples\\hpc_examples.py").read())
        try:
            file_path = "/home/upb/stekiner/Bachelorarbeit/FEMMT/FEM_Magnetics_Toolbox/femmt/examples/hpc_examples.py"

            with open(file_path, "r") as file:
                exec(file.read())

        except Exception as e:
            traceback.print_exc()

        Combination = Combination + 1
