"""Test different winding options."""
# python libraries
import os
from typing import Dict, List


import femmt as fmt
inductor_combinations = [
    {
        "Name": "Single RoundSolid Hexagonal",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RoundSolid,
        "ConductorArrangement": fmt.ConductorArrangement.Hexagonal
    },
    {
        "Name": "Single RoundSolid Square",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RoundSolid,
        "ConductorArrangement": fmt.ConductorArrangement.Square
    },
    {
        "Name": "Single RoundSolid SquareFullWidth",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RoundSolid,
        "ConductorArrangement": fmt.ConductorArrangement.SquareFullWidth
    },
    {
        "Name": "Single RectangularSolid Full",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RectangularSolid,
        "WindingScheme": fmt.WindingScheme.Full
    },
    {
        "Name": "Single RectangularSolid FoilHorizontal",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RectangularSolid,
        "WindingScheme": fmt.WindingScheme.FoilHorizontal
    },
    {
        "Name": "Single RectangularSolid FoilVertical (fixed thickness)",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RectangularSolid,
        "WindingScheme": fmt.WindingScheme.FoilVertical,
        "WrapParaType": fmt.WrapParaType.FixedThickness
    },
    {
        "Name": "Single RectangularSolid FoilVertical (interpolate)",
        "WindingType": fmt.WindingType.Single,
        "ConductorType": fmt.ConductorType.RectangularSolid,
        "WindingScheme": fmt.WindingScheme.FoilVertical,
        "WrapParaType": fmt.WrapParaType.Interpolate
    },
]

transformer_combinations = [
    {
        "Name": "Interleaved RoundSolid HorizontalAlternating",
        "WindingType": fmt.WindingType.TwoInterleaved,
        "ConductorType": fmt.ConductorType.RoundSolid,
        "WindingScheme": fmt.InterleavedWindingScheme.HorizontalAlternating,
        "ConductorArrangement": fmt.ConductorArrangement.Square
    },
    {
        "Name": "Interleaved RoundLitz VerticalStacked (square)",
        "WindingType": fmt.WindingType.TwoInterleaved,
        "ConductorType": fmt.ConductorType.RoundLitz,
        "WindingScheme": fmt.InterleavedWindingScheme.VerticalStacked,
        "ConductorArrangement": fmt.ConductorArrangement.Square
    },
    {
        "Name": "Interleaved RoundLitz VerticalStacked (hexagonal)",
        "WindingType": fmt.WindingType.TwoInterleaved,
        "ConductorType": fmt.ConductorType.RoundLitz,
        "WindingScheme": fmt.InterleavedWindingScheme.VerticalStacked,
        "ConductorArrangement": fmt.ConductorArrangement.Hexagonal
    },
]

def run_inductor_simulations(working_directory: str, combinations: List[Dict]):
    """
    Run the simulations to test several winding options for the inductor.

    :param working_directory: working directory
    :type working_directory: str
    :param combinations: combinations to simulate in a dictionary which are stored in a list
    :type combinations: List
    """
    not_working = []
    for combination in combinations:
        geo = fmt.MagneticComponent(fmt.ComponentType.Inductor, working_directory, True)

        core_db = fmt.core_database()["PQ 40/40"]
        core = fmt.Core(core_db["core_inner_diameter"], core_db["window_w"], core_db["window_h"],
                        "95_100") 
        geo.set_core(core)
        
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 50, 0.0005)
        geo.set_air_gaps(air_gaps)
            
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.001, 0.001)
        insulation.add_winding_insulations([0.0001, 0.0001], 0.001)
        geo.set_insulation(insulation)

        conductor = fmt.Conductor(0, fmt.Conductivity.Copper)

        conductor_arrangement = combination["ConductorArrangement"] if "ConductorArrangement" in combination else False

        if combination["ConductorType"] == fmt.ConductorType.RoundSolid:
            conductor.set_solid_round_conductor(0.0013, conductor_arrangement)
        elif combination["ConductorType"] == fmt.ConductorType.RoundLitz:
            conductor.set_litz_round_conductor(0.0013, 150, 100e-6, None, conductor_arrangement)
        else:
            conductor.set_rectangular_conductor(0.0013)

        winding_window = fmt.WindingWindow(core, insulation)
        complete = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
        winding_scheme = combination["WindingScheme"] if "WindingScheme" in combination else None
        wrap_para_type = combination["WrapParaType"] if "WrapParaType" in combination else None
        complete.set_winding(conductor, 15, winding_scheme, wrap_para_type)

        geo.set_winding_window(winding_window)

        name = combination["Name"]
        try:
            geo.create_model(freq=250000, pre_visualize_geometry=False, save_png=True)
        except Exception as e:
            print(e)
            not_working.append(name + ": " + str(e))

        image_path = os.path.join(working_directory, "mesh", "hybrid_color.png")
        os.rename(image_path, os.path.join(working_directory, "..", "images", f"{name}.png"))

    print(not_working)


def run_transformer_simulations(working_directory: str, combinations: List[Dict]):
    """
    Run the simulations to test several winding options for the transformer.

    :param working_directory: working directory
    :type working_directory: str
    :param combinations: combinations to simulate in a dictionary which are stored in a list
    :type combinations: List
    """
    not_working = []
    for combination in combinations:
        geo = fmt.MagneticComponent(fmt.ComponentType.Transformer, working_directory, True)

        core_db = fmt.core_database()["PQ 40/40"]
        core = fmt.Core(core_db["core_inner_diameter"], core_db["window_w"], core_db["window_h"],
                        "95_100") 
        geo.set_core(core)
        
        air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
        air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 50, 0.0005)
        geo.set_air_gaps(air_gaps)
            
        insulation = fmt.Insulation()
        insulation.add_core_insulations(0.001, 0.001, 0.001, 0.001)
        insulation.add_winding_insulations([0.0001, 0.0001], 0.001)
        geo.set_insulation(insulation)

        conductor1 = fmt.Conductor(0, fmt.Conductivity.Copper)
        conductor2 = fmt.Conductor(1, fmt.Conductivity.Copper)

        conductor_arrangement = combination["ConductorArrangement"] if "ConductorArrangement" in combination else False

        if combination["ConductorType"] == fmt.ConductorType.RoundSolid:
            conductor1.set_solid_round_conductor(0.0013, conductor_arrangement)
            conductor2.set_solid_round_conductor(0.0013, conductor_arrangement)
        elif combination["ConductorType"] == fmt.ConductorType.RoundLitz:
            conductor1.set_litz_round_conductor(0.0013, 150, 100e-6, None, conductor_arrangement)
            conductor2.set_litz_round_conductor(0.0013, 150, 100e-6, None, conductor_arrangement)
        else:
            conductor1.set_rectangular_conductor(0.0013)
            conductor2.set_rectangular_conductor(0.0013)

        winding_window = fmt.WindingWindow(core, insulation)
        if combination["WindingType"] == fmt.WindingType.TwoInterleaved:
            complete = winding_window.split_window(fmt.WindingWindowSplit.NoSplit)
            complete.set_interleaved_winding(conductor1, 10, conductor2, 10, combination["WindingScheme"], 0.0005)
        else:
            top, bot = winding_window.split_window(fmt.WindingWindowSplit.HorizontalSplit)
            winding_scheme = combination["WindingScheme"] if "WindingScheme" in combination else None
            wrap_para_type = combination["WrapParaType"] if "WrapParaType" in combination else None
            top.set_winding(conductor1, 5, winding_scheme, wrap_para_type)
            bot.set_winding(conductor2, 5, winding_scheme, wrap_para_type)

        geo.set_winding_window(winding_window)

        name = combination["Name"]
        try:
            geo.create_model(freq=250000, pre_visualize_geometry=True, save_png=True)
        except Exception as e:
            print(e)
            not_working.append(name + ": " + str(e))

        image_path = os.path.join(working_directory, "mesh", "hybrid_color.png")
        os.rename(image_path, os.path.join(working_directory, "..", "images", f"{name}.png"))

    print(not_working)


if __name__ == "__main__":
    working_directory = os.path.join(os.path.dirname(__file__), "winding_test")
    if not os.path.isdir(os.path.join(os.path.dirname(__file__), "images")):
        os.mkdir = os.path.join(os.path.dirname(__file__), "images")

    # run_inductor_simulations(working_directory, [inductor_combinations[2]])
    run_transformer_simulations(working_directory, transformer_combinations[0:3])
