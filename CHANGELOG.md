# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]


## [0.5.2] - 2024-04-30
### Added 
- log_material.json output for material logging information 
### Fixed
- Improved non-linear solver, as there were some one-shot results without any iteration
- Wrong displayed currents when using excitation_sweep() with complex currents
- Fixed incomplete python package in version 0.5.1

## [0.5.1] - 2024-02-06
### Fixed
- Fix documentation issues
- Fix material database dependency issues

## [0.5.0] - 2024-02-06
### Added
- various integration tests, unit tests for materialdatabase
- three winding transformer
- center-tapped transformer drawing schemes
- stacked transformer
- parallel connection of solid turns

### Changed
- API has lots of changes. Check out the examples and the documentation.

### Updated
- materialdatabase: material loading and interpolation of operation point

## [0.4.0] - 2022-12-19
### Added
- updated and improved syntax
- cost functions for core and wire material according to IEEE paper 'Component cost models for multi-objective optimizations of switched-mode power converters'
- connect femmt to the new pip package for the material database
- add optimization routine for automated design process
- Update GUI according to material database connection and optimization routine


## [0.3.0] - 2022-09-01
### Added
- Add color dictionaries for individual geometry visualization
- Added output json file for thermal simulation
- Add Parser to read and visualize the result.json-files
- Add first version of reluctance model including an example file
- Added dynamic mesh density algorithm for the winding window
- Simulation settings are now stored in the log file. A simulation can be started using given log-file.
- Added a new interface for femmt

### Fixed
- fix #13: improve reading the onelab filepath 
- fix #16: Wrong Mesh is simulated when changing the number of turns
- fix #17: Error with integrated_transformer mesh
- fix #19: Scale update in result plots
- fix #22: Fix bug for air-gap positions 0 and 100 percent

## [0.2.1] - 2022-04-28
### Updated
- possibility to assign fixed magnetic loss angle and conductivity in update_core function 
- new isolation scheme (looking from core to all windings): core_cond_isolation=[top, bottom, inner, outer] instead of core_cond_isolation=[prim2core, sec2core]
- gmsh 4.9.5 as minimum requirement

### Added
- example for foil winding inductor
- conductor material can be chosen from a small material database, used in update_conductors(), e.g. conductivity_sigma=["copper"]
- isolations for thermal simulation

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

[Unreleased]: https://github.com/upb-lea/transistordatabase/compare/0.5.1...HEAD
[0.5.1]: https://github.com/upb-lea/transistordatabase/compare/0.5.1...0.5.0
[0.5.0]: https://github.com/upb-lea/transistordatabase/compare/0.5.0...0.4.0
[0.4.0]: https://github.com/upb-lea/transistordatabase/compare/0.4.0...0.3.0
[0.3.0]: https://github.com/upb-lea/transistordatabase/compare/0.3.0...0.2.1
[0.2.1]: https://github.com/upb-lea/transistordatabase/compare/0.2.0...0.2.1
[0.2.0]: https://github.com/upb-lea/transistordatabase/compare/0.1.2...0.2.0
[0.1.2]: https://github.com/upb-lea/transistordatabase/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/upb-lea/transistordatabase/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/upb-lea/transistordatabase/compare/0.1.0...0.1.0


