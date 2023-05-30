import femmt as fmt
import os
import time
import shutil

def copy_electro_magnetic_necessary_files(src_folder, dest_folder):
    files = ["core_materials_temp.pro", "fields.pro", "ind_axi_python_controlled.pro", "Parameter.pro", "postquantities.pro", "solver.pro", "values.pro"]

    for file in files:
        from_path = os.path.join(src_folder, file)
        to_path = os.path.join(dest_folder, file)
        shutil.copy(from_path, to_path)

def create_example_model(working_directory, electro_magnetic_base_folder, strands_coefficients_folder_path):
    inductor_frequency = 270000
    
    # Copy important files from electro magnetic base folder
    electro_magnetic_folder = os.path.join(working_directory, "electro_magnetic")
    try:
        if not os.path.isdir(working_directory):
            os.mkdir(working_directory)
        if not os.path.isdir(electro_magnetic_folder):
            os.mkdir(electro_magnetic_folder)
        copy_electro_magnetic_necessary_files(electro_magnetic_base_folder, electro_magnetic_folder)
    except OSError as err:
        print(err)

    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, 
                                verbosity=fmt.Verbosity.ToFile, 
                                electro_magnetic_folder_path=electro_magnetic_folder,
                                strands_coefficients_folder_path=strands_coefficients_folder_path)
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
                    permittivity_measurement_setup="LEA_LK", silent=True)
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
    #winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6,fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)
    winding.parallel = False
    vww.set_winding(winding, 9, None)
    geo.set_winding_windows([winding_window])

    return geo

if __name__ == "__main__":
    electro_magnetic_base_folder = os.path.join(os.path.dirname(__file__), "..", "electro_magnetic")
    strands_coefficients_folder_path = os.path.join(electro_magnetic_base_folder, "Strands_Coefficients")
    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    benchmark_results_folder = os.path.join(example_results_folder, "benchmarks")
    parallel_folder = os.path.join(benchmark_results_folder, "parallel")

    if not os.path.exists(parallel_folder):
        os.mkdir(parallel_folder)

    inductor_frequency = 270000

    geos = []
    simulation_parameters = []
    working_directories = []
    for i in range(10):
        working_directory = os.path.join(parallel_folder, f"inductor_{i}")
        geos.append(create_example_model(working_directory, electro_magnetic_base_folder, strands_coefficients_folder_path))
        simulation_parameters.append({
            "freq": inductor_frequency,
            "current": [4.5]
        })

    start_time = time.time()
    fmt.hpc(5, geos, simulation_parameters)
    execution_time = time.time() - start_time

    print(f"Execution time: {execution_time}")

