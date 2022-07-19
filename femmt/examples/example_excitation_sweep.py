# minimal example for github clone
import femmt as fmt
import numpy as np

#component = "inductor"
component = "transformer"

# Create Object
if component == "inductor":

    # 0: choose frequencies, amplitude and phases to sweep
    frequencies = [100000, 200000]
    current_amplitudes = [[10], [4]]
    phases = [[0], [179]]

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor)

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]

    core = fmt.Core(core_w=core_db["core_w"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, None, 0.0005)
    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters: use solid wires
    winding = fmt.Winding(8, 0, fmt.Conductivity.Copper, fmt.WindingType.Primary, fmt.WindingScheme.Square)
    winding.set_litz_conductor(None, 100, 70e-6, 0.5)
    geo.set_windings([winding])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0001)
    geo.set_isolation(isolation)

    # 5. create the model
    geo.create_model(freq=100000, visualize_before=True, save_png=False)

    # 6. start simulation
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases, show_last=True)


if component == "transformer":
    # 0: choose frequencies, amplitude and phases to sweep
    frequencies = [100000, 200000]
    current_amplitudes = [[4, 14.5], [2, 6]]
    phases = [[0, 176], [0, 163]]

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer)

    # 2. set core parameters
    core = fmt.Core(window_h=0.0295, window_w=0.012, core_w=0.015,
                    non_linear=False, sigma=1, re_mu_rel=3200, phi_mu_deg=10)

    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 50, 0.0005)
    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters: use solid wires
    winding1 = fmt.Winding(21, 0, fmt.Conductivity.Copper, fmt.WindingType.Primary, fmt.WindingScheme.Square)
    winding1.set_litz_conductor(None, 100, 70e-6, 0.5)

    winding2 = fmt.Winding(0, 7, fmt.Conductivity.Copper, fmt.WindingType.Secondary, fmt.WindingScheme.Square)
    winding2.set_litz_conductor(None, 100, 70e-6, 0.5)

    geo.set_windings([winding1, winding2])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0002, 0.0002, 0.0005)
    geo.set_isolation(isolation)

    # 5. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, visualize_before=True)

    # 6. start simulation
    geo.excitation_sweep(frequency_list=frequencies, current_list_list=current_amplitudes, phi_deg_list_list=phases, show_last=True)