"""Planar matrix transformer with flexible split (3 columns × 2 rows) and interleaved windings."""
import femmt as fmt
import os
from datetime import datetime

Date = datetime.now().strftime("%Y%m%d-%H%M%S")

def interleaved_flexible_split_transformer(onelab_folder: str = None,
                                           show_visual_outputs: bool = True,
                                           is_test: bool = False):
    """Create a planar matrix transformer using a flexible split.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """

    def example_thermal_simulation(show_thermal_visual_outputs: bool = True, flag_insulation: bool = True):
        # Thermal simulation:
        # The losses calculated by the magnetics simulation can be used to calculate the heat distribution of the
        # given magnetic component. In order to use the thermal simulation, thermal conductivities for each material
        # can be entered as well as a boundary temperature which will be applied on the boundary of the
        # simulation (dirichlet boundary condition).

        # The case parameter sets the thermal conductivity for a case which will be set around the core.
        # This could model some case in which the transformer is placed in together with a set potting material.
        thermal_conductivity_dict = {
            "air": 0.0263,
            "case": {  # epoxy resign
                "top": 1.54,
                "top_right": 1.54,
                "right": 1.54,
                "bot_right": 1.54,
                "bot": 1.54
            },
            "core": 5,  # ferrite
            "winding": 400,  # copper
            "air_gaps": 180,  # aluminium nitride
            "insulation": 0.42 if flag_insulation else None  # polyethylene
        }

        # Here the case size can be determined
        case_gap_top = 0.002
        case_gap_right = 0.0025
        case_gap_bot = 0.002

        # Here the boundary temperatures can be set, currently it is set to 20°C (around 293°K).
        # This does not change the results of the simulation (at least when every boundary is set equally)
        # but will set the temperature offset.
        boundary_temperatures = {
            "value_boundary_top": 20,
            "value_boundary_top_right": 20,
            "value_boundary_right_top": 20,
            "value_boundary_right": 20,
            "value_boundary_right_bottom": 20,
            "value_boundary_bottom_right": 20,
            "value_boundary_bottom": 20
        }

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
        # When the losses file is already created and contains the losses for the current model, it is enough to
        # run geo.create_model in order for the thermal simulation to work (geo.single_simulation is not needed).
        # Obviously when the model is modified and the losses can be out of date and therefore the
        # geo.single_simulation needs to run again.
        geo.thermal_simulation(thermal_conductivity_dict, boundary_temperatures, boundary_flags, case_gap_top,
                               case_gap_right, case_gap_bot, show_thermal_visual_outputs,
                               color_scheme=fmt.colors_ba_jonas, colors_geometry=fmt.colors_geometry_ba_jonas,
                               flag_insulation=flag_insulation)

    # --- Setup working directories
    base_dir = os.path.dirname(__file__)
    results_dir = os.path.join(base_dir, "example_results")
    os.makedirs(results_dir, exist_ok=True)
    work_dir = os.path.join(results_dir, f"flexible_matrix_{Date}")
    os.makedirs(work_dir, exist_ok=True)

    # --- Initialize MagneticComponent
    geo = fmt.MagneticComponent(
        component_type=fmt.ComponentType.Transformer,
        working_directory=work_dir,
        is_gui=is_test
    )
    if onelab_folder:
        geo.file_data.onelab_folder_path = onelab_folder

    # --- 1) Core setup
    core_dims = fmt.dtos.SingleCoreDimensions(core_inner_diameter=5e-3,
                                                    window_w=6e-3,
                                                    window_h=2e-3,
                                                    core_h=8e-3)
    core_material = fmt.LinearComplexCoreMaterial(mu_r_abs=3100,
                                                  phi_mu_deg=12,
                                                  dc_conductivity=0.6,
                                                  eps_r_abs=0,
                                                  phi_eps_deg=0,
                                                  mdb_verbosity=fmt.Verbosity.Silent)
    core = fmt.Core(material=core_material,
                    core_type=fmt.CoreType.Single,
                    core_dimensions=core_dims,
                    detailed_core_model=False)
    geo.set_core(core)

    # --- 2) Air-gap setup
    gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(gaps)

    # --- 3) Insulation setup
    prepreg = 0.00023
    half_iso = prepreg / 2
    ins = fmt.Insulation(flag_insulation=True)
    ins.add_core_insulations(0.0005, 0.0005, 0.0005, 0.0005)
    ins.add_winding_insulations([[half_iso] * 2, [half_iso] * 2])
    geo.set_insulation(ins)

    # --- 4) Winding window & flexible split
    ww = fmt.WindingWindow(core, ins)
    h_cuts = [1/3, 2/3]  # split at 1/3 and 2/3 width → 3 columns
    v_cuts = [[0.5], [0.5], [0.5]]  # each column split at half height → 2 rows per col
    cells = ww.flexible_split(horizontal_split_factors=h_cuts,
                              vertical_split_factors=v_cuts)
    # cells is a flat list of 6 VirtualWindingWindow objects in order:
    # [TL, BL, TM, BM, TR, BR] if row-major per column
    # Empirically, use index mapping:
    top_left = cells[0]
    bottom_left = cells[1]
    top_mid = cells[2]
    bottom_mid = cells[3]
    top_right = cells[4]
    bottom_right = cells[5]

    # --- 5) Define conductors
    wrap_type = fmt.WrapParaType.FixedThickness
    p = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
    p.set_rectangular_conductor(thickness=33e-6)
    s = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
    s.set_rectangular_conductor(thickness=35e-6)

    # --- 6) Assign interleaved windings
    # Top row: P, S, P
    for region, cond in zip([top_left, top_mid, top_right], [p, s, p]):
        region.set_winding(cond, 1, fmt.WindingScheme.FoilHorizontal,
                           wrap_para_type=wrap_type,
                           foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)
    # Bottom row: S, P, S
    for region, cond in zip([bottom_left, bottom_mid, bottom_right], [s, p, s]):
        region.set_winding(cond, 1, fmt.WindingScheme.FoilHorizontal,
                           wrap_para_type=wrap_type,
                           foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)

    geo.set_winding_windows([ww])

    # --- 7) Simulation
    geo.create_model(freq=1e6, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(freq=1e6,
                          current=[7, 7],
                          phi_deg=[0, 180],
                          show_fem_simulation_results=show_visual_outputs)
    example_thermal_simulation(show_visual_outputs, flag_insulation=True)


if __name__ == "__main__":
    interleaved_flexible_split_transformer(show_visual_outputs=True)
