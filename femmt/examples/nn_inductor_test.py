"""Basic example to show how to simulate an inductor."""
import femmt as fmt
import materialdatabase as mdb
import os
import pandas as pd
from tabulate import tabulate

def inductor_simulation(onelab_folder: str = None, show_visual_outputs: bool = True, is_test: bool = False):
    """Run the example code for the inductor."""

    df = pd.read_pickle('test_df')
    print(df.shape)
    print(tabulate(df[df["Turns"]>12], headers='keys', tablefmt='pretty'))

    selected_rows = [83842]
    #selected_rows = [69, 428, 700]

    # df = df[(df["Total Core Eddy Losses_%-Error"] >= 100) & (df["Total Core Eddy Losses_%-Error"] >= 99)]
    print(df.shape)
    for index, row_idx in enumerate(selected_rows):
        row = df.loc[row_idx]

        frequency = row['Frequency']
        current = row['Current']
        turns = int(row['Turns'])
        temperature = row['Temperature']
        corediameter = row['Core Diameter']
        coreheight = row['Core Height']
        windowheight = row['Window Height']
        windowwidth = row['Window Width']
        airgaps = row['Air Gap Number']
        airgap_height = row['Air Gap Height']
        airgap_position = [row['Air Gap 1 Position'], row['Air Gap 2 Position'], row['Air Gap 3 Position']]
        coreinsulation_left = row['Core Insulation Left Core']
        conductor_radius = row['Conductor Radius']
        coductor_arrangement_hexagonal = row['Conductor_Arrangement_Hexagonal']
        coductor_arrangement_square = row['Conductor_Arrangement_Square']
        coductor_arrangement_squarefullwidth = row['Conductor_Arrangement_SquareFullWidth']

        print(f"Index: {index}")
        print(tabulate([row.to_dict()], headers="keys", tablefmt="grid"))
        print("\n")

        if coductor_arrangement_hexagonal == True:
            conductor_arrangement = fmt.ConductorArrangement.Hexagonal
        if coductor_arrangement_square == True:
            conductor_arrangement = fmt.ConductorArrangement.Square
        if coductor_arrangement_squarefullwidth == True:
            conductor_arrangement = fmt.ConductorArrangement.SquareFullWidth


        example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
        if not os.path.exists(example_results_folder):
            os.mkdir(example_results_folder)

        # Working directory can be set arbitrarily
        working_directory = os.path.join(example_results_folder, "inductor")
        if not os.path.exists(working_directory):
            os.mkdir(working_directory)

        # 1. chose simulation type
        geo = fmt.MagneticComponent(simulation_type=fmt.SimulationType.FreqDomain, component_type=fmt.ComponentType.Inductor, working_directory=working_directory,
                                    verbosity=fmt.Verbosity.ToConsole, is_gui=is_test)

        inductor_frequency = frequency

        # 2. set core parameters
        core_db = fmt.core_database()["PQ 40/40"]
        core_dimensions = fmt.dtos.SingleCoreDimensions(core_inner_diameter= corediameter,
                                                        window_w=windowwidth,
                                                        window_h=windowheight,
                                                        core_h=coreheight)

        core = fmt.Core(core_type=fmt.CoreType.Single,
                        core_dimensions=core_dimensions,
                        detailed_core_model=False,
                        material=mdb.Material.N95, temperature=temperature, frequency=inductor_frequency,
                        # permeability_datasource="manufacturer_datasheet",
                        permeability_datasource=fmt.MaterialDataSource.Measurement,
                        permeability_datatype=fmt.MeasurementDataType.ComplexPermeability,
                        permeability_measurement_setup=mdb.MeasurementSetup.LEA_LK,
                        permittivity_datasource=fmt.MaterialDataSource.Measurement,
                        permittivity_datatype=fmt.MeasurementDataType.ComplexPermittivity,
                        permittivity_measurement_setup=mdb.MeasurementSetup.LEA_LK, mdb_verbosity=fmt.Verbosity.Silent)

        # mu_rel=3000, phi_mu_deg=10,
        # sigma=0.5)
        geo.set_core(core)

        # 3. set air gap parameters
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        if airgaps in [1, 2, 3]:
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, airgap_height,
                                 airgap_position[0])
        # sets second air gap
        if airgaps in [2, 3]:
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, airgap_height,
                                 airgap_position[1])
        # sets third air gap
        if airgaps in [3]:
            air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, airgap_height,
                                 airgap_position[2])
        geo.set_air_gaps(air_gaps)

        # 4. set insulations
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, coreinsulation_left, 0.001)
        insulation.add_winding_insulations([[0.0005]])
        geo.set_insulation(insulation)

        # 5. create winding window and virtual winding windows (vww)
        winding_window = fmt.WindingWindow(core, insulation)
        vww = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)

        # 6. create conductor and set parameters: use solid wires
        winding = fmt.Conductor(0, fmt.Conductivity.Copper)
        winding.set_solid_round_conductor(conductor_radius=conductor_radius, conductor_arrangement = conductor_arrangement)
        winding.parallel = False


        # 7. add conductor to vww and add winding window to MagneticComponent
        vww.set_winding(winding, turns, None)
        geo.set_winding_windows([winding_window])

        # 8. create the model
        geo.create_model(freq=inductor_frequency, pre_visualize_geometry = True,save_png=False)

        # 6.a. start simulation
        geo.single_simulation(freq=inductor_frequency, current=[current],
                              plot_interpolation=False, show_fem_simulation_results=show_visual_outputs)




if __name__ == "__main__":
    inductor_simulation(show_visual_outputs=True)
