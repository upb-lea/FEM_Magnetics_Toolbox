import femmt as fmt
import os

component = "inductor_sweep"
# component = "transformer_sweep"

example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

# Create Object
if component == "inductor_sweep":

    # 0: choose frequencies, amplitude and phases to sweep
    frequencies = [100000, 200000]
    current_amplitudes = [[10], [4]]
    phases = [[0], [179]]

    # Working directory can be set arbitrarily
    working_directory = os.path.join(example_results_folder, "inductor_sweep")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory, silent=True)

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]

    core = fmt.Core(core_inner_diameter=core_db["core_inner_diameter"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                    material="N95", temperature=25, frequency=100000,
                    permeability_datasource=fmt.MaterialDataSource.Measurement,
                    permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                    permeability_measurement_setup="LEA_LK",
                    permittivity_datasource=fmt.MaterialDataSource.Measurement,
                    permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                    permittivity_measurement_setup="LEA_LK")

    # mu_rel=3000, phi_mu_deg=10,
    # sigma=0.5)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 10)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 90)
    geo.set_air_gaps(air_gaps)

    # 4. set insulations
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.004, 0.001)
    insulation.add_winding_insulations([0.0005])
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    # 6. create conductor and set parameters: use solid wires
    winding = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding.set_solid_round_conductor(conductor_radius=0.0013, conductor_arrangement=fmt.ConductorArrangement.Square)
    #winding.set_litz_round_conductor(conductor_radius=0.0013, number_strands=150, strand_radius=100e-6, 
    # fill_factor=None, conductor_arrangement=fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    vww.set_winding(winding, 9, None)
    geo.set_winding_window(winding_window)

    # 8. create the model
    geo.create_model(freq=100000, visualize_before=True, save_png=False)

    # 9. start simulation
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases, show_last=True)


if component == "transformer_sweep":
    # 0: choose frequencies, amplitude and phases to sweep
    frequencies = [100000, 200000]
    current_amplitudes = [[4, 14.5], [2, 6]]
    phases = [[0, 176], [0, 163]]

    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, "transformer_sweep")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer, working_directory=working_directory, silent=True)

    # 2. set core parameters
    core = fmt.Core(window_h=0.0295, window_w=0.012, core_inner_diameter=0.015,
                    mu_r_abs=3100, phi_mu_deg=12,
                    sigma=1.2)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 50, 0.0005)
    geo.set_air_gaps(air_gaps)

    # 4. set insulation
    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([0.0002, 0.0002], 0.0005)
    geo.set_insulation(insulation)

    # 5. create winding window and virtual winding windows (vww)
    winding_window = fmt.WindingWindow(core, insulation)
    left, right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit)

    # 6. create conductors and set parameters
    winding1 = fmt.Conductor(0, fmt.Conductivity.Copper)
    winding1.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    winding2 = fmt.Conductor(1, fmt.Conductivity.Copper)
    winding2.set_solid_round_conductor(0.0011, fmt.ConductorArrangement.Square)

    # 7. add conductor to vww and add winding window to MagneticComponent
    left.set_winding(winding1, 10, None)
    right.set_winding(winding2, 10, None)
    geo.set_winding_window(winding_window)

    # 8. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, visualize_before=True)

    # 9. start simulation
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases, show_last=True)