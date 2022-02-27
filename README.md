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

#### FEMMT release version (recommended)
This is the stable release version.
```
pip install femmt
```

#### FEMMT development version (for developers)
This is the latest development version with the latest features.
Note: You may need to install [git](https://git-scm.com/downloads).
```
cd /Documents/Folder/of/Interest   
git clone git@github.com:upb-lea/FEM_Magnetics_Toolbox.git
pip install -e .
```

### 2.4 Minimal example and first run
Run the example from here: [basic_example.py](/femmt/examples/basic_example.py).
FEMMT will ask you for the installation path of ONELAB during first use.


## 3. Examples
This toolbox is able to build a complete FEM simulation from simple Python code. The following figure shows the Python code on the left and the corresponding FEM simulation on the right.
![](https://github.com/upb-lea/FEM_Magnetics_Toolbox/blob/main/documentation/FEMMT_Screenshot.png?raw=true)

Code examples can be found in this [example file](/femmt/examples/basic_example.py). The examples contain among other things:
* Geometries: Coil, Transformer, Transformer with integrated stray path.
* Wire and stranded wire
* Air gaps
* Excitation with different frequencies

The simulation results can be found in the file /python-side-packages-path/femmt/femmt/results/result_log_electro_magnetic.json. In it you can find
* power loss in the core: hysteresis losses and eddy current losses
* Losses per winding and for each individual winding
* Self- and mutual inductances

## 4. GUI
There is a first preview for a GUI. Installing this is a bit cumbersome at first, but will be simplified in the future:
* Download the complete repository via `Code` -> `Download ZIP` and unpack it.
* install the development version of femmt as described above
* run python `downloads/path-to_femmt/femmt/gui/femmt_gui.py`
![](https://github.com/upb-lea/FEM_Magnetics_Toolbox/blob/main/documentation/femmt_gui_definition.png.png?raw=true)

## 5. Roadmap
Planned features in 2021
* work with pre-defined standard core structures
* work with pre-defined standard core materials
* automated 2D FEM simulation controlled with python

Planned features in 2022
* software stability and general improvements
* basic GUI implementation
* implement basics for thermal simulation

## 6. Bug Reports
Please use the issues report button within github to report bugs.

## 7. Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
For contributing, please refer to this [section](Contributing.md).

## 8. Changelog
Find the changelog [here](CHANGELOG.md)

## 9. License
[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

## 10. History and project status
This project was initially written in matlab using femm simulation tool. It became clear that the project was no longer a small project. The project should be completely rewritten, because many new complex levels have been added. To place the project in the open source world, the programming language python is used.      
