from femmt.femmt_enumerations import AirGapLegPosition, AirGapMethod
import femmt as fmt
import numpy as np
import os

def example_thermal_simulation():
    # Thermal simulation:
    # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the given magnetic component
    # In order to use the thermal simulation, thermal conductivities for each material can be entered as well as a boundary temperature
    # which will be applied on the boundary of the simulation (dirichlet boundary condition).
    
    # The case parameter sets the thermal conductivity for a case which will be set around the core.
    # This could model some case in which the transformer is placed in together with a set potting material.
    thermal_conductivity_dict = {
            "air": 0.0263,
            "case": { # (epoxy resign) | transformer oil
                "top": 0.122,
                "top_right": 0.122,
                "right": 0.122,
                "bot_right": 0.122,
                "bot": 0.122
            },
            "core": 5, # ferrite
            "winding": 400, # copper
            "air_gaps": 180, # aluminiumnitride
            "isolation": 0.42 # polyethylen
    }

    # Here the case size can be determined
    case_gap_top = 0.002
    case_gap_right = 0.0025
    case_gap_bot = 0.002

    # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
    # This does not change the results of the simulation (at least when every boundary is set equally) but will set the temperature offset.
    boundary_temperatures = {
        "value_boundary_top": 20,
        "value_boundary_top_right": 20,
        "value_boundary_right_top": 20,
        "value_boundary_right": 20,
        "value_boundary_right_bottom": 20,
        "value_boundary_bottom_right": 20,
        "value_boundary_bottom": 20
    }
    
    # In order to compare the femmt thermal simulation with a femm heat flow simulation the same boundary temperature should be applied.
    # Currently only one temperature can be applied which will be set on every boundary site.
    femm_boundary_temperature = 20

    # Here the boundary sides can be turned on (1) or off (0)
    # By turning off the flag a neumann boundary will be applied at this point with heat flux = 0
    boundary_flags = {
        "flag_boundary_top": 0,
        "flag_boundary_top_right": 0,
        "flag_boundary_right_top": 1,
        "flag_boundary_right": 1,
        "flag_boundary_right_bottom": 1,
        "flag_boundary_bottom_right": 1,
        "flag_boundary_bottom": 1
    }

    # In order for the thermal simulation to work an electro_magnetic simulation has to run before.
    # The em-simulation will create a file containing the losses.
    # When the losses file is already created and contains the losses for the current model, it is enough to run geo.create_model in
    # order for the thermal simulation to work (geo.single_simulation is not needed).
    # Obviously when the model is modified and the losses can be out of date and therefore the geo.single_simulation needs to run again.
    geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                           case_gap_right, case_gap_bot, True, color_scheme=fmt.colors_ba_jonas, colors_geometry=fmt.colors_geometry_ba_jonas)

    # Because the isolations inside of the winding window are not implemented in femm simulation.
    # The validation only works when the isolations for the FEMMT thermal simulation are turned off.
    geo.femm_thermal_validation(thermal_conductivity_dict, femm_boundary_temperature, case_gap_top, case_gap_right, case_gap_bot)

component = "inductor"
# component = "transformer-interleaved"
# component = "transformer"
# component = "integrated_transformer"

# Create Object
if component == "inductor":
    # Working directory can be set arbitrarily
    #working_directory = os.path.join(os.path.dirname(__file__), "working_directory")

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor) #working_directory=working_directory)

    # 2. set core parameters
    core_db = fmt.core_database()["PQ 40/40"]

    core = fmt.Core(core_w=core_db["core_w"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Manually, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.001, 0.0005)
    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters: use solid wires
    winding = fmt.Winding(8, 0, fmt.Conductivity.Copper, fmt.WindingType.Primary, fmt.WindingScheme.Square)
    winding.set_solid_conductor(0.0015)
    geo.set_windings([winding])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0001)
    geo.set_isolation(isolation)

    # 5. create the model
    geo.create_model(freq=100000, visualize_before=True, save_png=False)

    # 6. start simulation
    geo.single_simulation(freq=100000, current=[3], show_results=True)

    # 7. prepare and start thermal simulation
    example_thermal_simulation()

    # Excitation Sweep Example
    # Perform a sweep using more than one frequency
    # fs = [0, 10000, 30000, 60000, 100000, 150000]
    # amplitude_list = [[10], [2], [1], [0.5], [0.2], [0.1]]
    # phase_list = [[0], [10], [20], [30], [40], [50]]
    # geo.excitation_sweep(frequency_list=fs, current_list_list=amplitude_list, phi_deg_list_list=phase_list)

    # Reference simulation using FEMM
    # geo.femm_reference(freq=100000, current=[1], sigma_cu=58, sign=[1], non_visualize=0)

if component == "transformer-interleaved":
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
    winding1 = fmt.Winding(21, 0, fmt.Conductivity.Copper, fmt.WindingType.Interleaved, fmt.WindingScheme.Horizontal)
    winding1.set_solid_conductor(0.0011)

    winding2 = fmt.Winding(0, 7, fmt.Conductivity.Copper, fmt.WindingType.Interleaved, fmt.WindingScheme.Horizontal)
    winding2.set_solid_conductor(0.0011)

    geo.set_windings([winding1, winding2])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0002, 0.0002, 0.0005)
    geo.set_isolation(isolation)

    # 5. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, visualize_before=True)
    geo.single_simulation(freq=250000, current=[4, 12], phi_deg=[0, 180], show_results=True)


    # other simulation options:
    # ------------------------
    # read inductances
    # geo.get_inductances(I0=8, op_frequency=250000, skin_mesh_factor=0.5)

    # perform a reference simulation using FEMM
    # geo.femm_reference(freq=250000, current=[4, 12], sign=[1, -1], non_visualize=0)

    example_thermal_simulation()
    

if component == "transformer":
    # Example for a transformer with multiple virtual winding windows.

    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Transformer)

    # 2. set core parameters
    core = fmt.Core(window_h=0.0295, window_w=0.012, core_w=0.015,
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)
    geo.set_core(core)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 50, 0.0005)
    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters
    winding1 = fmt.Winding(10, 0, fmt.Conductivity.Copper, fmt.WindingType.Primary, fmt.WindingScheme.Square)
    winding1.set_solid_conductor(0.0011)

    winding2 = fmt.Winding(0, 10, fmt.Conductivity.Copper, fmt.WindingType.Secondary, fmt.WindingScheme.Square)
    winding2.set_litz_conductor(None, 600, 35.5e-6, 0.6)

    geo.set_windings([winding1, winding2])

    # 5. set isolation
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0002, 0.0002, 0.0005)
    geo.set_isolation(isolation)

    # 6. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, visualize_before=True)
    geo.single_simulation(freq=250000, current=[4.14723021, 14.58960019], phi_deg=[- 1.66257715/np.pi*180, 170])

if component == "integrated_transformer":
    # 1. chose simulation type
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.IntegratedTransformer)

    # 2. set core parameters
    core = fmt.Core(window_h=0.03, window_w=0.011, core_w=0.02,
                    mu_rel=3100, phi_mu_deg=12,
                    sigma=0.6)
    geo.set_core(core)

    # 2.1 set stray path parameters
    stray_path = fmt.StrayPath(start_index=0, radius=geo.core.core_w / 2 + geo.core.window_w - 0.001, width=None, midpoint=None)
    geo.set_stray_path(stray_path)

    # 3. set air gap parameters
    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 30, 0.001)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 40, 0.001)
    geo.set_air_gaps(air_gaps)

    # 4. set conductor parameters
    winding1 = fmt.Winding(1, 3, fmt.Conductivity.Copper, fmt.WindingType.Interleaved, fmt.WindingScheme.Horizontal)
    winding1.set_litz_conductor(None, 100, 70e-6, 0.5)

    winding2 = fmt.Winding(2, 6, fmt.Conductivity.Copper, fmt.WindingType.Interleaved, fmt.WindingScheme.Horizontal)
    winding2.set_litz_conductor(None, 100, 70e-6, 0.5)
    
    geo.set_windings([winding1, winding2])

    # 5. set isolations
    isolation = fmt.Isolation()
    isolation.add_core_isolations(0.001, 0.001, 0.002, 0.001)
    isolation.add_winding_isolations(0.0002, 0.0002, 0.0005)
    geo.set_isolation(isolation)

    # 6. start simulation with given frequency, currents and phases
    geo.create_model(freq=250000, visualize_before=True)
    geo.single_simulation(freq=250000, current=[8.0, 4.0], phi_deg=[0, 180])

    # other simulation options:
    # -------------------------
    # geo.get_inductances(I0=10, op_frequency=100000, skin_mesh_factor=0.5)
