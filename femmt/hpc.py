from multiprocessing import Pool
from femmt import MagneticComponent
from typing import List, Dict
import sys, os

def hpc_inner(parameters: Dict):
    model = parameters["model"]
    freq = parameters["freq"]
    current = parameters["current"]

    model.create_model(freq=freq)
    model.single_simulation(freq=freq, current=current)

def hpc(n_processes, models: List[MagneticComponent], simulation_parameters: List[Dict]):
    with Pool(processes=n_processes) as pool:
        parameters = []
        for index, (model, simulation_parameter) in enumerate(zip(models, simulation_parameters)):
            # Check simulation parameters
            if "freq" not in simulation_parameter:
                print(f"Missing simulation parameter {index}:freq. Simulation will be skipped.")
                continue
            if "current" not in simulation_parameter:
                print(f"Missing simulation parameter {index}:current. Simulation will be skipped.")
                continue

            parameters.append({
                "model": model,
                "freq": simulation_parameter["freq"],
                "current": simulation_parameter["current"]
            })
            
        pool.map(hpc_inner, parameters)

    print("Done")