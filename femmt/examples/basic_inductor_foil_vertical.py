import femmt as fmt
import os


def basic_example_inductor_foil_vertical(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    # Choose wrap para type
    wrap_para_type = fmt.WrapParaType.FixedThickness
    # wrap_para_type = fmt.WrapParaType.Interpolate

    example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
    if not os.path.exists(example_results_folder):
        os.mkdir(example_results_folder)

    # Example for a transformer with multiple virtual winding windows.
    working_directory = os.path.join(example_results_folder, "inductor_foil_vertical")
    if not os.path.exists(working_directory):
        os.mkdir(working_directory)

    # Choose wrap para type
    wrap_para_type = fmt.WrapParaType.FixedThickness

    # Set is_gui = True so FEMMt won't ask for the onelab path if no config is found.
    geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                verbosity=fmt.Verbosity.Silent, is_gui=True)

    core_db = fmt.core_database()["PQ 40/40"]
    core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter=core_db["core_inner_diameter"],
                                                    window_w=core_db["window_w"],
                                                    window_h=core_db["window_h"],
                                                    core_h=core_db["core_h"])

    core = fmt.Core(core_type=fmt.CoreType.Single, core_dimensions=core_dimensions,
                    mu_r_abs=3100, phi_mu_deg=12,
                    sigma=0.6, permeability_datasource=fmt.MaterialDataSource.Custom,
                    permittivity_datasource=fmt.MaterialDataSource.Custom)
    geo.set_core(core)

    air_gaps = fmt.AirGaps(fmt.AirGapMethod.Center, core)
    air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 0.0005)
    geo.set_air_gaps(air_gaps)

    insulation = fmt.Insulation()
    insulation.add_core_insulations(0.001, 0.001, 0.002, 0.001)
    insulation.add_winding_insulations([[0.0005]])
    geo.set_insulation(insulation)

    winding_window = fmt.WindingWindow(core, insulation)
    vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

    winding = fmt.Conductor(0, fmt.Conductivity.Copper, winding_material_temperature=25)
    winding.set_rectangular_conductor(thickness=1e-3)

    vww.set_winding(winding, 5, fmt.WindingScheme.FoilVertical, wrap_para_type)
    geo.set_winding_windows([winding_window])

    geo.create_model(freq=100000, pre_visualize_geometry=show_visual_outputs, save_png=False)

    geo.single_simulation(freq=100000, current=[3], show_fem_simulation_results=show_visual_outputs)


if __name__ == "__main__":
    basic_example_inductor_foil_vertical(show_visual_outputs=True)