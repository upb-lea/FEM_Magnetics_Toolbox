# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.2] - 2021-08-08
### Updated
- updated strand approximation
### Added
- add option for dedicated stray path
- add complex permeability and permitivity for core materials for Core Loss estimation
- add iGSE and GSE for Core Loss estimation

## [0.1.1] - 2021-08-11
### Updated
- updated inductance calculations
- code clean up

### Fixed
- fix #2: config.json was not read correct
- fix #3: Install pyfemm on windows machines in case of not installed pyfemm

## [0.1.0] - 2021-07-28
### Added
#### Structure
- add README.md
- add CHANGELOG.md
- add femmt/__init__.py

#### Essentials
- add femmt/FEMMT.py
- add femmt/functions.py
- add femmt/ind_axi_python_controlled.pro
- add femmt/solver.pro
- add femmt/BH.pro

#### Examples
- add femmt/FEMMT_geometric.py
- add femmt/basic_example.py

#### Additional/Experimental Code
- add femmt/pandas_json.py
- add femmt/femm_test.py
- add femmt/SimComparison.py
- add femmt/SolidComp.py
- add femmt/CompRes.py

[Unreleased]: https://github.com/upb-lea/FEM_Magnetics_Toolbox/compare/0.1.1...HEAD
[0.1.2]: https://github.com/upb-lea/FEM_Magnetics_Toolbox/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/upb-lea/FEM_Magnetics_Toolbox/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/upb-lea/FEM_Magnetics_Toolbox/compare/0.1.0...0.1.0



