# Python standard libraries
from multiprocessing import Pool
from typing import List, Dict, Callable
import os
import shutil

# Third parry libraries
from femmt import MagneticComponent
from femmt.data import FileData

def _copy_electro_magnetic_necessary_files(src_folder, dest_folder):
    files = ["core_materials_temp.pro", "fields.pro", "ind_axi_python_controlled.pro", "Parameter.pro", "postquantities.pro", "solver.pro", "values.pro"]

    for file in files:
        from_path = os.path.join(src_folder, file)
        to_path = os.path.join(dest_folder, file)
        shutil.copy(from_path, to_path)

def _hpc(parameters: Dict):
    model = parameters["model"]
    freq = parameters["simulation_parameters"]["freq"]
    current = parameters["simulation_parameters"]["current"]

    model.create_model(freq=freq, pre_visualize_geometry=False, save_png=False)
    model.single_simulation(freq=freq, current=current, plot_interpolation=False, show_fem_simulation_results=False)

def run(n_processes: int, models: List[MagneticComponent], simulation_parameters: List[Dict], working_directory: str, custom_hpc: Callable = None):
    """Executes the given models on the given number of parallel processes. Typically this number shouldn't be higher than the number of cores of the processor.

    :param n_processes: Number of parallel processes. If this is equal to None the number returned by os.cpu_count() is used.
    :type n_processes: int
    :param models: List of MagneticComponents which shall be simulated in parallel.
    :type models: List[MagneticComponent]
    :param simulation_parameters: List of dictionaries containing the parameters for the parallel simulation. The Nth item corresponds to the Nth MagneticComponent in models.
                                    For the default hpc_function this dictionary needs the frequency and the current (which is needed for the single_simulation function).
                                    Since the dictionary is given to the hpc_function directly, if a custom_hpc function is used the needed parameters can be added through this dict.
    :type simulation_parameters: List[Dict]
    :param working_directory: The directory which will store the model data and results data for every parallel simulation.
    :type working_directory: str
    :param custom_hpc: If set to None the default hpc will be used (create_model() and single_simulation() are executed). If a custom_hpc is set this function
                        will be called for the parallel execution and the funtion parameters is the simulation_parameter dict taken from the simulation_parameters list., defaults to None
    :type custom_hpc: Callable, optional
    """
    electro_magnetic_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "electro_magnetic")
    strands_coefficients_folder = os.path.join(electro_magnetic_folder, "Strands_Coefficients")
    
    if not os.path.isdir(working_directory):
        os.mkdir(working_directory)

    for index, model in enumerate(models):
        # Setup necessary files and directories
        model_name = model.simulation_name if model.simulation_name is not None else f"model_{index}"
        model_working_directory = os.path.join(working_directory, model_name)
        model_electro_magnetic_directory = os.path.join(model_working_directory, "electro_magnetic")
        if not os.path.isdir(model_working_directory):
            os.mkdir(model_working_directory)
        if not os.path.isdir(model_electro_magnetic_directory):
            os.mkdir(model_electro_magnetic_directory)

        _copy_electro_magnetic_necessary_files(electro_magnetic_folder, model_electro_magnetic_directory)

        # Update directories for each model
        model.file_data.update_paths(model_working_directory, model_electro_magnetic_directory, strands_coefficients_folder)

    # Create pool of workers and apply _hpc to it
    with Pool(processes=n_processes) as pool:
        parameters = []
        for index, (model, simulation_parameter) in enumerate(zip(models, simulation_parameters)):
            # Check simulation parameters
            if custom_hpc is None:
                if "freq" not in simulation_parameter:
                    print(f"Missing simulation parameter {index}:freq. Simulation will be skipped.")
                    continue
                if "current" not in simulation_parameter:
                    print(f"Missing simulation parameter {index}:current. Simulation will be skipped.")
                    continue

            parameters.append({
                "model": model,
                "simulation_parameters": simulation_parameter
            })
            
        if custom_hpc is None:
            pool.map(_hpc, parameters)
        else:
            pool.map(custom_hpc, parameters)