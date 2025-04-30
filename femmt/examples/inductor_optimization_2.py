import optuna
import femmt as fmt
import numpy as np
import os
import json

# inductance_range = (0.0002, 0.0009)  # H

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
    working_directory = os.path.join(example_results_folder, "inductor_parasitic_optimization")
    trial_folder = f"trial_r{int(conductor_radius * 1e6)}_n{n_turns}_s{scheme}"
    trial_working_dir = os.path.join(working_directory, trial_folder)
    os.makedirs(trial_working_dir, exist_ok=True)

    core_db = fmt.core_database()[core_name]
    bobbin_db = fmt.bobbin_database()[core_name]

    core_dimensions = fmt.dtos.SingleCoreDimensions(
        core_inner_diameter=core_db["core_inner_diameter"],
        window_w=core_db["window_w"],
        window_h=core_db["window_h"],
        core_h=core_db["core_h"]
    )
    bobbin_dimensions = fmt.dtos.BobbinDimensions(
        bobbin_inner_diameter=bobbin_db["bobbin_inner_diameter"],
        bobbin_window_w=bobbin_db["bobbin_window_w"],
        bobbin_window_h=bobbin_db["bobbin_window_h"],
        bobbin_h=bobbin_db["bobbin_h"]
    )

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
    winding.set_solid_round_conductor(conductor_radius, fmt.ConductorArrangement.Square)
    vww.set_winding(winding, int(n_turns), None, fmt.Align.ToEdges, placing_strategy=fmt.ConductorDistribution(scheme))
    geo.set_winding_windows([winding_window])

    try:
        geo.create_model(freq=270000, pre_visualize_geometry=False, save_png=False)
        geo.single_simulation(freq=270000, current=[1], show_fem_simulation_results=False)

        inductance = read_flux_over_current(trial_working_dir)
        # if not (inductance_range[0] <= inductance <= inductance_range[1]):
        #     print(f"⚠️ Skipped: Inductance {inductance:.4e} H out of range.")
        #     return float("-inf"), float("inf")  # prefer lowest possible inductance, high capacitance

        capacitance = geo.get_inductor_stray_capacitance(show_fem_simulation_results=False)
        return inductance, capacitance

    except Exception as e:
        print(f"❌ Simulation failed: {e}")
        return float("-inf"), float("inf")

def objective(trial):
    core_name = trial.suggest_categorical("core_name", ["PQ 40/40"])
    conductor_radius = 0.225e-3
    air_gap_height = 0.001
    n_turns = 84
    scheme = trial.suggest_int("winding_scheme", 1, 8)
    return run_simulation(core_name, conductor_radius, n_turns, air_gap_height, scheme)

if __name__ == "__main__":
    study = optuna.create_study(directions=["maximize", "minimize"])
    study.optimize(objective, n_trials=1000)

    import optuna.visualization
    fig = optuna.visualization.plot_pareto_front(
        study,
        target_names=["Inductance", "Capacitance"]
    )
    fig.show()
