# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Updated

### Added

### Fixed
- fix #15: Secondary to Core isolation thickness not working


## [0.2.0] - 2022-02-14
### Updated
- updated the winding generation
- updated the femm reference model
- updated project structure
- updated meshing to hybrid mesh
- updated class structure
### Added
- add file Analytical_Core_Data.py
- add example files DAB_Input_Data.py, DAB_trafo_optimization.py
- add file mu_imag.pro
- add folder femmt/thermal
- add result_log_electro_magnetic
- add horizontal interleaved winding scheme
- add method write_log() to femmt.py
- add thermal simulation with onelab
- add femm heat flow validation with femm
### Fixed
- fix #5: changed typo to L_h_conc = self.M**2 / self.L_22
- fix #11: rename femmt.py to femmt_classes.py due to package problems

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

[Unreleased]: https://github.com/upb-lea/transistordatabase/compare/0.2.0...HEAD
[0.2.0]: https://github.com/upb-lea/transistordatabase/compare/0.1.2...0.2.0
[0.1.2]: https://github.com/upb-lea/transistordatabase/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/upb-lea/transistordatabase/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/upb-lea/transistordatabase/compare/0.1.0...0.1.0


