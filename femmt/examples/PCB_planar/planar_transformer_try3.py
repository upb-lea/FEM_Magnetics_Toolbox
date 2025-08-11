"""Create a planar transformer with 2-4 turns using the virtual winding window."""
import femmt as fmt
import os
from datetime import datetime
import logging

Date = datetime.now().strftime("%Y%m%d-%H%M%S")
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def transformer_3x2(onelab_folder: str = None,
                    show_visual_outputs: bool = True,
                    is_test: bool = False):
    """
    Planar matrix transformer with: P | S | P.

    :param onelab_folder: onelab folder path
    :type onelab_folder: str
    :param show_visual_outputs: True to show visual outputs (simulation results)
    :type show_visual_outputs: bool
    :param is_test: True for pytest usage. Defaults to False.
    :type is_test: bool
    """
    # 1) Setup output folder
    base = os.path.dirname(__file__)
    results = os.path.join(base, "example_results")
    os.makedirs(results, exist_ok=True)
    work = os.path.join(results, f"3x2_{Date}")
    os.makedirs(work, exist_ok=True)

    # 2) Initialize magnetic component
    geo = fmt.MagneticComponent(
        component_type=fmt.ComponentType.Transformer,
        working_directory=work,
        is_gui=is_test
    )
    if onelab_folder:
        geo.file_data.onelab_folder_path = onelab_folder

    # 3) Core definition
    dims = fmt.dtos.SingleCoreDimensions(
        core_inner_diameter=5e-3,
        window_w=6e-3,
        window_h=2e-3,
        core_h=8e-3
    )
    core = fmt.Core(
        core_type=fmt.CoreType.Single,
        core_dimensions=dims,
        mu_r_abs=3100,
        phi_mu_deg=12,
        sigma=0.6,
        permeability_datasource=fmt.MaterialDataSource.Custom,
        permittivity_datasource=fmt.MaterialDataSource.Custom
    )
    geo.set_core(core)

    # 4) Air gap in center leg for tunable leakage Lk
    ag = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    ag.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(ag)

    # 5) Insulation setup
    prep = 0.00023
    ins = fmt.Insulation(flag_insulation=True)
    ins.add_core_insulations(0.0005, 0.0005, 0.0005, 0.0005)
    ins.add_winding_insulations([[prep/2, prep/2], [prep/2, prep/2]])
    # Kapton material is added between every layer of turns
    layer_insulation = fmt.insulation_materials_database()["film_insulation"]["Kapton"]
    ins.add_insulation_between_layers(thickness=2.8e-4, dielectric_constant=layer_insulation["dielectric_constant"])
    geo.set_insulation(ins)

    # 6) Split into 3 columns × 2 rows
    ww = fmt.WindingWindow(core, ins)
    cells = ww.flexible_split(
        horizontal_split_factors=[1/3, 2/3],      # cuts at 33% and 67% → 3 columns
        vertical_split_factors=[[0.5], [0.5], [0.5]]  # each column splits at 50% → 2 rows
    )
    # flexible_split returns flat list length=6:
    # [col0-row0, col0-row1,  col1-row0, col1-row1,  col2-row0, col2-row1]
    c0r0, c0r1, c1r0, c1r1, c2r0, c2r1 = cells

    # 7) Conductors
    wrap = fmt.WrapParaType.FixedThickness
    P = fmt.Conductor(0, fmt.Conductivity.Copper)
    P.set_rectangular_conductor(thickness=35e-6)
    S = fmt.Conductor(1, fmt.Conductivity.Copper)
    S.set_rectangular_conductor(thickness=35e-6)

    gap = 0.0002  # 0.2 mm clearance

    # 8) Assign each row: [P, S, P]
    # Top row: c0r0, c1r0, c2r0 → P, S, P (all routed downward)
    for region, cond in zip([c0r0, c1r0, c2r0], [P, S, P]):
        region.set_winding(
            cond, turns=1,
            winding_scheme=fmt.WindingScheme.FoilHorizontal,
            wrap_para_type=wrap,
            foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward
        )

    # Bottom row: c0r1, c1r1, c2r1 → P, S, P (all routed upward to clear)
    for region, cond in zip([c0r1, c1r1, c2r1], [P, S, P]):
        region.set_winding(
            cond, turns=1,
            winding_scheme=fmt.WindingScheme.FoilHorizontal,
            wrap_para_type=wrap,
            foil_horizontal_placing_strategy=fmt.FoilHorizontalDistribution.VerticalDownward
        )

    geo.set_winding_windows([ww])

    # 9) Run EM simulation
    geo.create_model(freq=1e6, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(
        freq=1e6,
        current=[7, 7],
        phi_deg=[0, 180],
        show_fem_simulation_results=show_visual_outputs)


if __name__ == "__main__":
    transformer_3x2(show_visual_outputs=True)
