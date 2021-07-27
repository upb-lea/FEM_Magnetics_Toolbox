# FEM_Magnetics_Toolbox
Python toolbox to generate preconfigured figures for FEM simulation tools

Functionality examples
 * work with pre-defined standard core structures [not implemented yet]
 * use python to perform parametersweeps, e.g. perform several automated simultations of different air gap sizes
 * read the results automatet with python from the FEM simulation tool



## Installation
```
cd /Documents/Folder/of/Interest   
git clone git@github.com:upb-lea/FEM_Magnetics_Toolbox.git
```


## Basics of Usage
Date: 27.07.2021
How to use the FEM Magnetics Toolbox:
* Import the class with "from FEMMT import MagneticComponent"
* Create instance with for example "geo = MagneticComponent(component_type="transformer")"
				or "geo = MagneticComponent(component_type="inductor")"
				...
* Define/Update the Core geometry with "geo.update_core(core_type="EI", window_h=0.03)"
* Add/Update the various air gaps with "geo.update_air_gaps(method="percent", n_air_gaps=2, position_tag=[0, 0], air_gap_h=[0.001, 0.001], air_gap_position=[20, 80])"
				    or "geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.002])"
				    ...
* Define/Update conductor(s) with "geo.update_conductors(n_turns=[8, 12], conductor_type=["solid", "solid"], conductor_radix=[0.0015, 0.001])
* Start a single simulation with "geo.single_simulation(freq=100000, current=[5, 10], phi=[0, 0], skin_mesh_factor=accuracy)"

Mninimal example for a single simulation with displayed result in ONELAB:
```
from FEMMT import MagneticComponent

# Create Object
geo = MagneticComponent(component_type="inductor")

# Update Geometry
geo.update_core(core_type="EI", window_h=0.03)
geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.001])
geo.update_conductors(n_turns=[22], conductor_type=["solid"], conductor_radix=[0.0015])

# Perform a single simulation
geo.single_simulation(freq=1000000, current=[10])
```

## Roadmap
Planned features in 2021
* work with pre-defined standard core structures
* work with pre-defined standard core materials
* automated 2D FEM simulation controlled with python

## Bug Reports
Please use the issues report button within github to report bugs.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
For contributing, please refer to this [section](Contributing.md).

For drawings (e.g. in readme-files), we recomment to use the program [Inkscape](https://inkscape.org/). It is open source software and runs on Linux, Mac and Windows. If you want to draw electirc circuits, we recommend this library on [github](https://github.com/upb-lea/Inkscape_electric_Symbols).

## Authors and acknowledgement
Current developers
 * Till Piepenbrock
 * Jan Wiegard
 * Dennis Kirchner

Developers in the past


## Changelog
Find the changelog [here](CHANGELOG.md)

## License
[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## History and project status
This project was initially written in matlab using femm simulation tool. It became clear that the project was no longer a small project. The project should be completely rewritten, because many new complex levels have been added. To place the project in the open source world, the programming language python is used.      
