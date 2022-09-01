import femmt as fmt
import os

# Working directory can be set arbitrarily
working_directory = os.path.join(os.path.dirname(__file__), "inductor")
if not os.path.exists(working_directory):
    os.mkdir(working_directory)

geo = fmt.MagneticComponent(component_type=fmt.ComponentType.Inductor, working_directory=working_directory)

core_db = fmt.core_database()["PQ 40/40"]

core = fmt.Core(core_w=core_db["core_w"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                material="95_100")

geo.set_core(core)

air_gaps = fmt.AirGaps(fmt.AirGapMethod.Percent, core)
air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 10, 0.0005)
air_gaps.add_air_gap(fmt.AirGapLegPosition.CenterLeg, 90, 0.0005)
geo.set_air_gaps(air_gaps)

isolation = fmt.Isolation()
isolation.add_core_isolations(0.001, 0.001, 0.004, 0.001)
isolation.add_winding_isolations(0.0005)
geo.set_isolation(isolation)

winding_window = fmt.WindingWindow(core, isolation, 0.5, 0.5)
top_left, top_right, bot_left, bot_right = winding_window.virtual_winding_windows

left = winding_window.combine_vww(top_left, bot_left)
right = winding_window.combine_vww(top_right, bot_right)

complete = winding_window.combine_vww(left, right)

conductor1 = fmt.Conductor(0, fmt.Conductivity.Copper)
conductor1.set_solid_round_conductor(0.0013, fmt.ConductorArrangement.Square)

complete.set_winding(conductor1, 9, None)
geo.set_virtual_winding_windows([complete])

geo.create_model(freq=100000, visualize_before=True, save_png=False)
geo.single_simulation(freq=100000, current=[4.5], show_results=True)
