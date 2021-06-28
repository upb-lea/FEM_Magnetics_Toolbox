# FEM_Magnetics_Toolbox
Python toolbox to generate preconfigured figures for FEM simulation tools

Functionality examples
 * work with pre-defined standard core structures
 * use python to perform parametersweeps, e.g. perform several automated simultations of different air gap sizes
 * read the results automatet with python from the FEM simulation tool



## Installation
```
cd /Documents/Folder/of/Interest   
git clone git@github.com:upb-lea/FEM_Magnetics_Toolbox.git
```


## Usage
Date: 04.05.2021
How to use the FEM Magnetics Toolbox:
* Import the class with "from FEMMT import MagneticComponent"
* Createte instance with for example "geo = MagneticComponent()"
* Call a method with for example "geo.single_simulation() or geo.freq_sweep_simulation(start=0, end=250000, steps=6)"
Mninimal example for a single simulation with displayed result in ONELAB:
```
from FEMMT import MagneticComponent
import numpy as np

# Create Object
geo = MagneticComponent()
# Perform a single simulation
geo.single_simulation(freq=100000,current=1)
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
