import matplotlib.pyplot as plt
import femmt as fmt
import os

if not os.path.exists(os.path.join(os.path.dirname(__file__), "sweep_examples")):
    os.mkdir(os.path.join(os.path.dirname(__file__), "sweep_examples"))

# Change between different example sweeps
sweep = "air_gap_height"

if sweep == "air_gap_height":
    # In this sweep an inductor with variable air gap height is simulated
    air_gap_heights = [0.0020, 0.0010, 0.0005, 0.000025, 0.0000125]

    working_directories = []

    for i, height in enumerate(air_gap_heights):
        # In order to save the results for every simulation the working directory is changed
        directory = os.path.join(os.path.dirname(__file__), 'sweep_examples', f"air_gap_{i}")
        if not directory:
            os.mkdir(directory)

        working_directories.append(directory)

        geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=directory, silent=True)

        core_db = fmt.core_database()["PQ 40/40"]
        core = fmt.Core(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                        material="95_100")
        geo.set_core(core)

        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, height, 50)
        geo.set_air_gaps(air_gaps)

        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
        insulation.add_winding_insulations([0.0005])
        geo.set_insulation(insulation)

        winding_window = fmt.WindingWindow(core, insulation)
        complete = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        conductor = fmt.Conductor(0, fmt.Conductivity.Copper)
        conductor.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6, 
                                            fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)
        complete.set_winding(conductor, 9, None)
        geo.set_winding_window(winding_window)

        geo.create_model(freq=100000, visualize_before=False, save_png=False)
        geo.single_simulation(freq=100000, current=[4.5], show_results=False)

    # After the simulations the sweep can be analyzed
    # This could be done using the FEMMTLogParser:
    logs = fmt.FEMMTLogParser.get_log_files_from_working_directories(working_directories)
    log_parser = fmt.FEMMTLogParser(logs)

    # In this case the self inductivity of winding1 will be analyzed
    inductivities = []
    for name, data in log_parser.data.items():
        inductivities.append(data.sweeps[0].windings[0].self_inductance)

    plt.plot(air_gap_heights, inductivities, "ro")
    plt.title("Air gap height sweep")
    plt.xlabel("air gap height")
    plt.ylabel("self inductance")
    plt.grid()
    plt.show()