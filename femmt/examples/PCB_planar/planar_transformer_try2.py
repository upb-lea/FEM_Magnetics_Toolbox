"""Planar matrix transformer with flexible split (2 columns × 3 rows) and custom winding positions."""
import femmt as fmt
import os
from datetime import datetime

Date = datetime.now().strftime("%Y%m%d-%H%M%S")

def customized_split_transformer(onelab_folder: str = None,
                                 show_visual_outputs: bool = True,
                                 is_test: bool = False):
    """Create a planar matrix transformer with a 2×3 grid of virtual winding windows.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    # --- Prepare directories
    base_dir = os.path.dirname(__file__)
    results_dir = os.path.join(base_dir, "example_results")
    os.makedirs(results_dir, exist_ok=True)
    work_dir = os.path.join(results_dir, f"custom_matrix_{Date}")
    os.makedirs(work_dir, exist_ok=True)

    # --- Initialize MagneticComponent
    geo = fmt.MagneticComponent(
        component_type=fmt.ComponentType.Transformer,
        working_directory=work_dir,
        is_gui=is_test
    )
    if onelab_folder:
        geo.file_data.onelab_folder_path = onelab_folder

    # --- Core definition
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

    # --- Air gap setup
    ag = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    ag.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(ag)

    # --- Insulation
    prepreg = 0.00023
    half_iso = prepreg/2
    ins = fmt.Insulation(flag_insulation=True)
    ins.add_core_insulations(0.0005, 0.0005, 0.0005, 0.0005)
    ins.add_winding_insulations([[half_iso]*2, [half_iso]*2])
    geo.set_insulation(ins)

    # --- Winding window & flexible split
    ww = fmt.WindingWindow(core, ins)
    # two columns (cut at 0.5) × three rows (cuts at 1/3, 2/3)
    h_cuts = [0.5]
    v_cuts = [[1/3, 2/3], [1/3, 2/3]]
    cells = ww.flexible_split(
        horizontal_split_factors=h_cuts,
        vertical_split_factors=v_cuts
    )
    # cells flat list order: [TL, ML, BL, TR, MR, BR]
    top_left, mid_left, bot_left, top_right, mid_right, bot_right = cells

    # --- Conductors
    wrap_type = fmt.WrapParaType.FixedThickness
    p = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
    p.set_rectangular_conductor(thickness=35e-6)
    s = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
    s.set_rectangular_conductor(thickness=35e-6)

    # --- Assign windings per request
    # Primary on top-left & top-right
    for region in [top_left, top_right]:
        region.set_winding(p, 1, fmt.WindingScheme.FoilHorizontal,
                           wrap_para_type=wrap_type,
                           foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)
    # Secondary on mid-left & mid-right
    for region in [mid_left, mid_right]:
        region.set_winding(s, 1, fmt.WindingScheme.FoilHorizontal,
                           wrap_para_type=wrap_type,
                           foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)
    # Primary on bot-left & bot-right
    for region in [bot_left, bot_right]:
        region.set_winding(p, 1, fmt.WindingScheme.FoilHorizontal,
                           wrap_para_type=wrap_type,
                           foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward)

    geo.set_winding_windows([ww])

    # --- Simulation
    geo.create_model(freq=1e6, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(freq=1e6,
                          current=[7, 7],
                          phi_deg=[0, 180],
                          show_fem_simulation_results=show_visual_outputs)


if __name__ == "__main__":
    customized_split_transformer(show_visual_outputs=True)
