"""Create a planar transformer with four turns using the virtual winding window."""
import femmt as fmt
import os
from datetime import datetime
import logging

Date = datetime.now().strftime("%Y%m%d-%H%M%S")
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def transformer_2x4(onelab_folder: str = None,
                    show_visual_outputs: bool = True,
                    is_test: bool = False):
    """Planar matrix transformer with P, S, P, S.

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
    work = os.path.join(results, f"2x4_{Date}")
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

    # 4) Air gap in center leg
    ag = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
    ag.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005, 50)
    geo.set_air_gaps(ag)

    # 5) Insulation setup
    prep = 0.00023
    ins = fmt.Insulation(flag_insulation=True)
    ins.add_core_insulations(0.0005, 0.0005, 0.0005, 0.0005)
    ins.add_winding_insulations([[prep/2, prep/2], [prep/2, prep/2]])
    # layer_ins = fmt.insulation_materials_database()["film_insulation"]["Kapton"]
    # ins.add_insulation_between_layers(
    #     thickness=1.5e-4,
    #     dielectric_constant=layer_ins["dielectric_constant"]
    # )
    geo.set_insulation(ins)

    # 6) Split into 2 columns × 4 rows
    ww = fmt.WindingWindow(core, ins)
    cells = ww.flexible_split(
        horizontal_split_factors=[0.25, 0.5, 0.75],  # 4 rows
        vertical_split_factors=[
            [0.5],  # each row split once → 2 columns
            [0.5],
            [0.5],
            [0.5],
        ]
    )
    # unpack 8 regions in row‐major order:
    r0c0, r0c1, r1c0, r1c1, r2c0, r2c1, r3c0, r3c1 = cells

    # 7) Conductors
    wrap = fmt.WrapParaType.FixedThickness
    P = fmt.Conductor(0, fmt.ConductorMaterial.Copper)
    P.set_rectangular_conductor(thickness=70e-6)
    S = fmt.Conductor(1, fmt.ConductorMaterial.Copper)
    S.set_rectangular_conductor(thickness=70e-6)

    # 8) Assign P→S→P→S in each column
    for col in (0, 1):
        regs = [cells[col + 2*row] for row in range(4)]
        seq = [P, S, P, S]
        dirs = [
            fmt.FoilHorizontalDistribution.VerticalDownward,
            fmt.FoilHorizontalDistribution.VerticalDownward,
            fmt.FoilHorizontalDistribution.VerticalDownward,
            fmt.FoilHorizontalDistribution.VerticalDownward
        ]
        for region, cond, direction in zip(regs, seq, dirs):
            region.set_winding(
                cond, turns=1,
                winding_scheme=fmt.WindingScheme.FoilHorizontal,
                wrap_para_type=wrap,
                foil_horizontal_placing_strategy=direction
            )

    geo.set_winding_windows([ww])

    # 9) Run EM simulation
    geo.create_model(freq=1e6, pre_visualize_geometry=show_visual_outputs)
    geo.single_simulation(
        freq=1e6,
        current=[1, 1],
        phi_deg=[0, 180],
        show_fem_simulation_results=show_visual_outputs)


if __name__ == "__main__":
    transformer_2x4(show_visual_outputs=True)
