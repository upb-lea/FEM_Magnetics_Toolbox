"""
Basic example to show how to simulate an inductor.

After starting the program, the geometry dimensions are displayed. Verify this geometry, close the window, to continue the simulation.
After a short time, B-Field and winding losses simulation results are shown. Winding losses are shown as a colormap.
In the core, the magnitude B-Field in Tesla is shown. With the gmsh window, one can move the picture in the 3D way (not recommended).
If you close this window, the thermal simulation will be continued, if programmed. If true, the thermal heat distribution will be displayed.
To continue with the next simulation (or end the program), you need to close this window. All results are written to the result
folder .../femmt/examples/example_results/simulation_file_name/results/log_electro_magnetic.json. and .../results_thermal.json.
"""
import femmt as fmt
import os
import logging

# configure logging to show femmt terminal output
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

example_results_folder = os.path.join(os.path.dirname(__file__), "example_results")
if not os.path.exists(example_results_folder):
    os.mkdir(example_results_folder)

# Working directory can be set arbitrarily
working_directory = os.path.join(example_results_folder, os.path.splitext(os.path.basename(__file__))[0])
if not os.path.exists(working_directory):
    os.mkdir(working_directory)

file_path = "C:/Users/tpiepe/Repositories/FEM_Magnetics_Toolbox/femmt/examples/example_results/aaa_new_mdb/results/log_electro_magnetic.json"
file_name = "log_electro_magnetic.json"
file_path_dict = {file_name: file_path}
# After the simulations the sweep can be analyzed
# This could be done using the FEMMTLogParser:
log_parser = fmt.FEMMTLogParser(file_path_dict)

frequency = 0
current = 0
for _, data in log_parser.data.items():
    frequency = data.sweeps[0].frequency
    current = data.sweeps[0].windings[0].current.real


geo = fmt.MagneticComponent.decode_settings_from_log(file_path, working_directory)

geo.create_model(freq=frequency, pre_visualize_geometry=False, save_png=False)

geo.single_simulation(freq=frequency, current=[current], show_fem_simulation_results=False)
