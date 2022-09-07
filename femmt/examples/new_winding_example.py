from femmt import MagneticComponent as fmc
from femmt import Functions as ff
from femmt import Model as model
from femmt.Enumerations import *
import os

# Working directory can be set arbitrarily
working_directory = os.path.join(os.path.dirname(__file__), "inductor")
if not os.path.exists(working_directory):
    os.mkdir(working_directory)

geo = fmc.MagneticComponent(component_type=ComponentType.Inductor, working_directory=working_directory)

core_db = ff.core_database()["PQ 40/40"]

core = model.Core(core_w=core_db["core_w"], window_w=core_db["window_w"], window_h=core_db["window_h"],
                material="95_100")

geo.set_core(core)

air_gaps = model.AirGaps(AirGapMethod.Percent, core)
air_gaps.add_air_gap(AirGapLegPosition.CenterLeg, 10, 0.0005)
air_gaps.add_air_gap(AirGapLegPosition.CenterLeg, 90, 0.0005)
geo.set_air_gaps(air_gaps)

isolation = model.Isolation()
isolation.add_core_isolations(0.001, 0.001, 0.004, 0.001)
isolation.add_winding_isolations(0.0005)
geo.set_isolation(isolation)

winding_window = model.WindingWindow(core, isolation)
# top_left, top_right, bot_left, bot_right = winding_window.split_window(fmt.WindingWindowSplit.HorizontalAndVerticalSplit, 0.5, 0.5)
complete = winding_window.split_window(WindingWindowSplit.NoSplit)

conductor1 = model.Conductor(0, Conductivity.Copper)
conductor1.set_solid_round_conductor(0.0013, ConductorArrangement.Square)
complete.set_winding(conductor1, 9, None)

geo.set_winding_window(winding_window)

geo.create_model(freq=100000, visualize_before=True, save_png=False)
geo.single_simulation(freq=100000, current=[4.5], show_results=True)
