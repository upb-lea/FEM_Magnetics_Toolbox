# This file is used for the LEA Master Project Group - Speed up of FEMMT by using parallel processing
# It contains an example simulation for which the simulation time and the needed memory size shall be measured.
import femmt as fmt
import os
import time
import json
import numpy as np

class SingleBenchmark:
    
    # Benchmark times
    setup_time: float                           # Time to setup the FEMMT objects 
    high_level_geo_gen_time: float              # Time to create the model in gmsh
    generate_hybrid_mesh_time : float           # Time to generate the hybrid mesh
    generate_electro_magnetic_mesh_time: float  # Time to generate the electro magnetic mesh
    prepare_simulation_time: float              # Time to prepare the simulation (create GetDP files and set simulation options)
    real_simulation_time: float                 # Time to run the simulation in GetDP
    logging_time: float                         # Time for FEMMT to write the results in a log

    create_model_time: float                    # high_level_geo_gen_time + generate_hybrid_mesh_time
    simulation_time: float                      # generate_electro_magnetic_mesh_time + prepare_simulation_time + real_simulation_time + logging_time

    execution_time_: float                      # all measurements summed up

    def __init__(self, benchmark_dict: dict):
        self.update_benchmark(benchmark_dict)

    def update_benchmark(self, benchmark_dict):
        self.setup_time = benchmark_dict["setup_time"]
        self.high_level_geo_gen_time = benchmark_dict["high_level_geo_gen_time"]
        self.generate_hybrid_mesh_time = benchmark_dict["generate_hybrid_mesh_time"]
        self.generate_electro_magnetic_mesh_time = benchmark_dict["generate_electro_magnetic_mesh_time"]
        self.prepare_simulation_time = benchmark_dict["prepare_simulation_time"]
        self.real_simulation_time = benchmark_dict["real_simulation_time"]
        self.logging_time = benchmark_dict["logging_time"]

        self.create_model_time = self.high_level_geo_gen_time + self.generate_hybrid_mesh_time
        self.simulation_time = self.generate_electro_magnetic_mesh_time + self.prepare_simulation_time + self.real_simulation_time + self.logging_time

        self.execution_time = self.setup_time + self.create_model_time + self.simulation_time

    def to_dict(self):
        return {
            "setup_time": self.setup_time,
            "create_model_time": self.create_model_time,
            "simulation_time": self.simulation_time,
            "high_level_geo_gen_time": self.simulation_time,
            "generate_hybrid_mesh_time": self.generate_hybrid_mesh_time,
            "generate_electro_magnetic_mesh_time": self.generate_electro_magnetic_mesh_time,
            "prepare_simulation_time": self.prepare_simulation_time,
            "real_simulation_time": self.real_simulation_time,
            "logging_time": self.logging_time,
            "execution_time": self.execution_time
        }

class Benchmark:
    benchmarks: SingleBenchmark

    def __init__(self, benchmarks):
        self.benchmarks = benchmarks

    def to_json(self, file_path):
        output = {
            "averages": {
                "setup_time": np.mean([x.setup_time for x in self.benchmarks]),
                "create_model_time": np.mean([x.create_model_time for x in self.benchmarks]),
                "simulation_time": np.mean([x.simulation_time for x in self.benchmarks]),
                "high_level_geo_gen_time": np.mean([x.high_level_geo_gen_time for x in self.benchmarks]),
                "generate_hybrid_mesh_time": np.mean([x.generate_hybrid_mesh_time for x in self.benchmarks]),
                "generate_electro_magnetic_mesh_time": np.mean([x.generate_electro_magnetic_mesh_time for x in self.benchmarks]),
                "prepare_simulation_time": np.mean([x.prepare_simulation_time for x in self.benchmarks]),
                "real_simulation_time": np.mean([x.real_simulation_time for x in self.benchmarks]),
                "logging_time": np.mean([x.logging_time for x in self.benchmarks]),
                "execution_time": np.mean([x.execution_time for x in self.benchmarks]),
            },
            "variances": {
                "setup_time": np.var([x.setup_time for x in self.benchmarks]),
                "create_model_time": np.var([x.create_model_time for x in self.benchmarks]),
                "simulation_time": np.var([x.simulation_time for x in self.benchmarks]),
                "high_level_geo_gen_time": np.var([x.high_level_geo_gen_time for x in self.benchmarks]),
                "generate_hybrid_mesh_time": np.var([x.generate_hybrid_mesh_time for x in self.benchmarks]),
                "generate_electro_magnetic_mesh_time": np.var([x.generate_electro_magnetic_mesh_time for x in self.benchmarks]),
                "prepare_simulation_time": np.var([x.prepare_simulation_time for x in self.benchmarks]),
                "real_simulation_time": np.var([x.real_simulation_time for x in self.benchmarks]),
                "logging_time": np.var([x.logging_time for x in self.benchmarks]),
                "execution_time": np.var([x.execution_time for x in self.benchmarks]),
            },
            "benchmarks": [x.to_dict() for x in self.benchmarks] 
        }

        with open(file_path, "w") as fd:
            json.dump(output, fd, indent=2)

def benchmark_simulation(working_directory):

    if not os.path.exists(working_directory):
        os.mkdir(working_directory)
    inductor_frequency = 270000

    start_time = time.time()
    
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, silent=True)
    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"])
    core = fmt.Core(core_type=fmt.CoreType.Single,
                    core_dimensions=core_dimensions,
                    material="N49", temperature=45, frequency=inductor_frequency,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup="LEA_LK",
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup="LEA_LK")
    geo.set_core(core)
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(air_gaps)
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False
    vww.set_winding(winding, 9, None)
    geo.set_winding_windows([winding_window])

    setup_time = time.time() - start_time

    high_level_geo_gen_time, generate_hybrid_mesh_time = geo.create_model(freq=inductor_frequency, pre_visualize_geometry=False, save_png=False, benchmark=True)

    generate_electro_magnetic_mesh_time, prepare_simulation_time, real_simulation_time, logging_time = geo.single_simulation(freq=inductor_frequency, current=[4.5],
                          plot_interpolation=False, show_fem_simulation_results=False, benchmark=True)

    return {
        "setup_time": setup_time,
        "high_level_geo_gen_time": high_level_geo_gen_time,
        "generate_hybrid_mesh_time": generate_hybrid_mesh_time,
        "generate_electro_magnetic_mesh_time": generate_electro_magnetic_mesh_time,
        "prepare_simulation_time": prepare_simulation_time,
        "real_simulation_time": real_simulation_time,
        "logging_time": logging_time 
    }
    
def single_benchmark(benchmark_results_folder):
    benchmark_results_file = os.path.join(benchmark_results_folder, "benchmark_results_single.json")
    benchmarks = []

    for i in range(30):
        benchmark = benchmark_simulation(os.path.join(benchmark_results_folder, "inductor"))
        benchmarks.append(SingleBenchmark(benchmark))

    benchmark = Benchmark(benchmarks)
    benchmark.to_json(benchmark_results_file)

if __name__ == "__main__":
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    benchmark_results_folder = os.path.join(example_results_folder, "benchmarks")

    if not os.path.exists(benchmark_results_folder):
        os.mkdir(benchmark_results_folder)

    single_benchmark(benchmark_results_folder)