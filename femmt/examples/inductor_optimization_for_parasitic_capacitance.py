import optuna
import femmt as fmt
import numpy as np
import os
import json


# def read_total_losses(results_path):
#     """
#     Reads total_losses from log_electro_magnetic.json inside results/values
#     """
#     try:
#         json_path = os.path.join(results_path, "results", "log_electro_magnetic.json")
#         with open(json_path, "r") as f:
#             data = json.load(f)
#         return data["total_losses"]["total_losses"]
#     except Exception as e:
#         print(f"❌ Failed to read total_losses: {e}")
#         return float("inf")
#
# def read_flux_over_current(results_path):
#     try:
#         json_path = os.path.join(results_path, "results", "log_electro_magnetic.json")
#         with open(json_path, "r") as f:
#             data = json.load(f)
#         return data["single_sweeps"][0]["winding1"]["flux_over_current"][0]
#     except Exception as e:
#         print(f"❌ Failed to read flux_over_current: {e}")
#         return None
#
# def run_simulation(conductor_radius, n_turns):
#     """
#     Run inductor simulation and return loss and capacitance as objectives.
#     Each simulation is saved in a unique subfolder.
#     """
#     example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
#     if not os.path.exists(example_results_folder):
#         os.mkdir(example_results_folder)
#
#     # Create subfolder specifically for inductor parasitic optimization
#     working_directory = os.path.join(example_results_folder, "inductor_parasitic_optimization")
#     if not os.path.exists(working_directory):
#         os.mkdir(working_directory)
#     trial_folder = f"trial_r{int(conductor_radius * 1e6)}_n{n_turns}"
#     trial_working_dir = os.path.join(working_directory, trial_folder)
#     os.makedirs(trial_working_dir, exist_ok=True)
#
#     # Core setup
#     core_db = fmt.core_database()["PQ 40/40"]
#     core_dimensions = fmt.dtos.SingleCoreDimensions(
#         core_inner_diameter=core_db["core_inner_diameter"],
#         window_w=core_db["window_w"],
#         window_h=core_db["window_h"],
#         core_h=core_db["core_h"]
#     )
#
#     core = fmt.Core(
#         core_type=fmt.CoreType.Single,
#         core_dimensions=core_dimensions,
#         material=fmt.Material.N95,
#         temperature=45,
#         frequency=270000,
#         permeability_datasource=fmt.MaterialDataSource.Measurement,
#         permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
#         permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
#         permittivity_datasource=fmt.MaterialDataSource.Measurement,
#         permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
#         permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK,
#         mdb_verbosity=fmt.Verbosity.Silent
#     )
#
#     geo = fmt.MagneticComponent(
#         simulation_type=fmt.SimulationType.FreqDomain,
#         component_type=fmt.ComponentType.Inductor,
#         working_directory=working_directory,
#         verbosity=fmt.Verbosity.Silent
#     )
#
#     geo.set_core(core)
#
#     # Air gap
#     air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
#     air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 50)
#     geo.set_air_gaps(air_gaps)
#
#     # Insulation
#     insulation = fmt.Insulation(flag_insulation=True)
#     insulation.add_core_insulations(1.55e-3, 1.55e-3, 0.9e-3, 1.5e-4)
#     insulation.add_winding_insulations([[0.025e-3]])
#     insulation.add_conductor_air_conductor_insulation([
#         [1e-3] * 6,
#         [1e-3, 1e-3]
#     ])
#     insulation.add_kapton_insulation(add_kapton=True, thickness=0.5e-3)
#     geo.set_insulation(insulation)
#
#     # Winding
#     winding_window = fmt.WindingWindow(core, insulation)
#     vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
#
#     winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=45)
#     winding.set_solid_round_conductor(conductor_radius=conductor_radius, conductor_arrangement=fmt.ConductorArrangement.Square)
#     vww.set_winding(winding, int(n_turns), None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution.VerticalUpward_HorizontalRightward)
#     geo.set_winding_windows([winding_window])
#
#     try:
#         geo.create_model(freq=270000, pre_visualize_geometry=False, save_png=False)
#         geo.single_simulation(freq=270000, current=[1], show_fem_simulation_results=False)
#
#         total_loss = read_total_losses(working_directory)
#         capacitance = geo.get_inductor_capacitance(show_fem_simulation_results=False)
#
#         return total_loss, capacitance
#
#     except Exception as e:
#         print(f"❌ Simulation failed: {e}")
#         return float("inf"), float("inf")
#
#
# def objective(trial):
#     conductor_radius = trial.suggest_float("conductor_radius", 0.1e-3, 0.5e-3)
#     # n_turns = trial.suggest_int("n_turns", 10, 80)
#     # Suggest n_turns from 20 evenly spaced values between 10 and 80
#     # turns_options = np.linspace(10, 80, 20, dtype=int).tolist()
#     # n_turns = trial.suggest_categorical("n_turns", turns_options)
#     n_turns = 5
#
#
#     loss, capacitance = run_simulation(conductor_radius, n_turns)
#     return loss, capacitance
#
#
# if __name__ == "__main__":
#     study = optuna.create_study(directions=["minimize", "minimize"])
#     study.optimize(objective, n_trials=5)
#
#     # Visualize Pareto front
#     import optuna.visualization
#     fig = optuna.visualization.plot_pareto_front(study)
#     fig.show()

# Target inductance value (flux_over_current[0])
# TARGET_INDUCTANCE = 0.0004818280353049755
# INDUCTANCE_TOLERANCE = 1e-6  # Allowable deviation
# inductance_range = (0.0003, 0.0005)
inductance_range = (0.0002, 0.0009)


def read_total_losses(results_path):
    try:
        json_path = os.path.join(results_path, "results", "log_electro_magnetic.json")
        with open(json_path, "r") as f:
            data = json.load(f)
        return data["total_losses"]["total_losses"]
    except Exception as e:
        print(f"❌ Failed to read total_losses: {e}")
        return float("inf")

def read_flux_over_current(results_path):
    try:
        json_path = os.path.join(results_path, "results", "log_electro_magnetic.json")
        with open(json_path, "r") as f:
            data = json.load(f)
        return data["single_sweeps"][0]["winding1"]["flux_over_current"][0]
    except Exception as e:
        print(f"❌ Failed to read flux_over_current: {e}")
        return None

def run_simulation(core_name, conductor_radius, n_turns, air_gap_height, scheme):
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    os.makedirs(example_results_folder, exist_ok=True)

    working_directory = os.path.join(example_results_folder, "inductor_parasitic_optimization")
    os.makedirs(working_directory, exist_ok=True)

    trial_folder = f"trial_r{int(conductor_radius * 1e6)}_n{n_turns}"
    trial_working_dir = os.path.join(working_directory, trial_folder)
    os.makedirs(trial_working_dir, exist_ok=True)

    # Core setup
    # core_db = fmt.core_database()["PQ 40/40"]
    core_db = fmt.core_database()[core_name]
    core_dimensions = fmt.dtos.SingleCoreDimensions(
        core_inner_diameter=core_db["core_inner_diameter"],
        window_w=core_db["window_w"],
        window_h=core_db["window_h"],
        core_h=core_db["core_h"]
    )
    bobbin_db = fmt.bobbin_database()[core_name]
    bobbin_dimensions = fmt.dtos.BobbinDimensions(bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
                                                  bobbin_window_w=bobbin_db["bobbin_window_w"],
                                                  bobbin_window_h=bobbin_db["bobbin_window_h"],
                                                  bobbin_h=bobbin_db["bobbin_h"])

    core = fmt.Core(
        core_type=fmt.CoreType.Single,
        core_dimensions=core_dimensions,
        bobbin_dimensions=bobbin_dimensions,
        material=fmt.Material.N95,
        temperature=45,
        frequency=270000,
        permeability_datasource=fmt.MaterialDataSource.Measurement,
        permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
        permeability_measurement_setup=fmt.MeasurementSetup.LEA_LK,
        permittivity_datasource=fmt.MaterialDataSource.Measurement,
        permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
        permittivity_measurement_setup=fmt.MeasurementSetup.LEA_LK,
        mdb_verbosity=fmt.Verbosity.Silent
    )

    geo = fmt.MagneticComponent(
        simulation_type=fmt.SimulationType.FreqDomain,
        component_type=fmt.ComponentType.Inductor,
        working_directory=trial_working_dir,
        verbosity=fmt.Verbosity.Silent
    )
    geo.set_core(core)

    # air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    # air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 50)
    # geo.set_air_gaps(air_gaps)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, air_gap_height, 50)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation(flag_insulation=True)
    insulation.add_core_insulations(1.55e-3, 1.55e-3, 0.9e-3, 1.5e-4)
    insulation.add_winding_insulations([[0.025e-3]])
    insulation.add_conductor_air_conductor_insulation([
        [1.0837e-4] * 6,
        [1e-3, 1e-3]
    ])
    insulation.add_kapton_insulation(add_kapton=True, thickness=0.07e-3)
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=45)
    winding.set_solid_round_conductor(conductor_radius=conductor_radius, conductor_arrangement=fmt.ConductorArrangement.Square)
    vww.set_winding(winding, int(n_turns), None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution(scheme))
    geo.set_winding_windows([winding_window])

    try:
        geo.create_model(freq=270000, pre_visualize_geometry=False, save_png=False)
        geo.single_simulation(freq=270000, current=[1], show_fem_simulation_results=False)

        inductance = read_flux_over_current(trial_working_dir)
        # if inductance is None or abs(inductance - TARGET_INDUCTANCE) > INDUCTANCE_TOLERANCE:
        #     print(f"⚠️ Trial skipped: inductance {inductance:.4e} not within tolerance of target {TARGET_INDUCTANCE:.4e}")
        #     return float("inf"), float("inf")
        # Filter by inductance range
        if not (inductance_range[0] <= inductance <= inductance_range[1]):
            print(f"⚠️ Skipped: Inductance {inductance:.4e} H out of range.")
            return float("inf"), float("inf")
        print(f"inductance is:", {inductance} )
        total_loss = read_total_losses(trial_working_dir)
        # capacitance = geo.get_inductor_capacitance(show_fem_simulation_results=False)
        capacitance = geo.get_inductor_stray_capacitance(show_fem_simulation_results=False, show_visual_outputs=False)

        return total_loss, capacitance

    except Exception as e:
        print(f"❌ Simulation failed: {e}")
        return float("inf"), float("inf")

def objective(trial):
    cores = [
        #"PQ 20/20",
        #"PQ 26/25",
        #"PQ 32/30",
        "PQ 40/40",
        #"PQ 50/50"
    ]
    core_name = trial.suggest_categorical("core_name", cores)
    # conductor_radius = trial.suggest_float("conductor_radius", 0.225e-3, 0.5e-3)
    conductor_radius = 0.225e-3
    # air_gap_height = trial.suggest_float("air_gap_height", 0.0005, 0.001)
    # air_gap_height = [0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.0010]
    # air_gap_height = trial.suggest_categorical(
    #     "air_gap_height", [0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.0010]
    # )
    air_gap_height = 0.001
    n_turns = 84
    scheme = trial.suggest_int("winding_scheme", 1, 8)
    # n_turns = 42
    # n_turns = trial.suggest_int("n_turns", 10, 80)
    loss, capacitance = run_simulation(core_name, conductor_radius, n_turns, air_gap_height, scheme)
    return loss, capacitance

if __name__ == "__main__":
    study = optuna.create_study(directions=["minimize", "minimize"])
    study.optimize(objective, n_trials=10)

    import optuna.visualization
    # fig = optuna.visualization.plot_pareto_front(study)
    # fig.show()
    fig = optuna.visualization.plot_pareto_front(
        study,
        target_names=["Losses", "Capacitance"]
    )
    fig.show()

