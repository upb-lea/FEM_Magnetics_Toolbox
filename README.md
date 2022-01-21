# FEM Magnetics Toolbox (FEMMT)
Python toolbox to generate preconfigured figures for FEM simulation tools in power electronics

Functionality examples
 * work with pre-defined standard core structures [not implemented yet]
 * use python to perform parametersweeps, e.g. perform several automated simulations of different air gap sizes
 * read the results automated with python from the FEM simulation tool

__Note: Alpha Version!__

## 1. Detailed Documentation
Can be found [here](https://upb-lea.github.io/FEM_Magnetics_Toolbox/main/intro.html).

## 2. Installation

### 2.1 ONELAB installation
* Go to https://onelab.info/
* Download the Desktop Version for your OS (Windows, Linux or macOS)
* Unpack the software and remember the file path. This will be needed later when installing FEMMT.

### 2.2 FEMM installation [for Windows User only]
* Go to https://www.femm.info/wiki/Download
* Install FEMM as described
* FEMM can be used as an alternative FEM solver for 2D simulations

### 2.3 install FEMMT
Chose to install the development version of FEMMT or the release version.

#### Installation the latest FEMMT development version (for developers)
Note: You may need to install [git](https://git-scm.com/downloads).
```
cd /Documents/Folder/of/Interest   
git clone git@github.com:upb-lea/FEM_Magnetics_Toolbox.git
pip install -e .
```

### 2.4 Minimal example and first run
Run the example from here: [basic_example.py](/femmt/Examples/basic_example.py).
FEMMT will ask you for the installation path of ONELAB during first use.


## 3. Basics of Usage
![](https://github.com/upb-lea/FEM_Magnetics_Toolbox/blob/main/documentation/Transformer_Screenshot.png?raw=true)
Date: 2021-07-27
How to use the FEM Magnetics Toolbox:
* Import the class with "from FEMMT import MagneticComponent"
* Create instance with for example "geo = MagneticComponent(component_type="transformer")"
				or "geo = MagneticComponent(component_type="inductor")"
				...
* Define/Update the Core geometry with "geo.update_core(core_type="EI", window_h=0.03)"
* Add/Update the various air gaps with "geo.update_air_gaps(method="percent", n_air_gaps=2, position_tag=[0, 0], air_gap_h=[0.001, 0.001], air_gap_position=[20, 80])"
				    or "geo.update_air_gaps(method="center", n_air_gaps=1, air_gap_h=[0.002])"
				    ...
* Define/Update conductor(s) for a transformer with "geo.update_conductors(n_turns=[8, 12], conductor_type=["solid", "solid"], conductor_radix=[0.0015, 0.001])
			    or for an inductor with "geo.update_conductors(n_turns=[[14]], conductor_type=["solid"], conductor_radix=[0.0015], winding=["primary"], scheme=["square"], core_cond_isolation=[0.0005], cond_cond_isolation=[0.0001])
* Start a single simulation with "geo.single_simulation(freq=100000, current=[5, 10], phi=[0, 0], skin_mesh_factor=accuracy)"

See examples in [basic_example.py](/femmt/Examples/basic_example.py)

## 4. Roadmap
Planned features in 2021
* work with pre-defined standard core structures
* work with pre-defined standard core materials
* automated 2D FEM simulation controlled with python

Planned features in 2022
* software stability and general improvements
* basic GUI implementation
* implement basics for thermal simulation

## 5. Bug Reports
Please use the issues report button within github to report bugs.

## 6. Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
For contributing, please refer to this [section](Contributing.md).

For drawings (e.g. in readme-files), we recomment to use the program [Inkscape](https://inkscape.org/). It is open source software and runs on Linux, Mac and Windows. If you want to draw electirc circuits, we recommend this library on [github](https://github.com/upb-lea/Inkscape_electric_Symbols).

## 7. Changelog
Find the changelog [here](CHANGELOG.md)

## 8. License
[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## 9. History and project status
This project was initially written in matlab using femm simulation tool. It became clear that the project was no longer a small project. The project should be completely rewritten, because many new complex levels have been added. To place the project in the open source world, the programming language python is used.      
