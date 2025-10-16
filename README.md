# Commands

- Git submodules

`git submodule update --remote --init`

- CMake configure step

`cmake --preset msvc-ninja-configure`

- CMake build step

`cmake --build --preset msvc-ninja-build`

- CMake configure + build

`cmake --workflow --preset msvc-ninja-all`

- Set your MSVC and Qt path in CMakePresets.json before attempting to configure or build
- To compile outside of IDE use x64 Native Tools Command Prompt for VS

# pagmo4scgms

- Intel TBB and Boost libraries must be installed on the host system and detected by CMake (they're pagmo dependencies)
- Compiled pagmo4scgms needs to have pagmo.dll in its directory, pagmo.dll then needs Boost and TBB DLLs to be either in the env path variable or in the same directory
- Additionally pagmo also depends on Eigen3 (without it some components of Pagmo won't work), Eigen3 is already in the /deps/ directory, however in order for CMake to detect it the Eigen3Configure file needs to exist. This file is created by running Eigen3's CMake configure step.
- "namespace ppr", "SPO" and duplicated nlopt code were removed from the original pagmo4scgms library (see commit details)
- TODO:
    - Consider compiling TBB and Boost in /deps/ directory (possibly very long download and compile time?)
        - Somehow utilize existing convention of passing deps by CMake command line params? 
    - Solve issues with pagmo4scgms.dll's dependencies - prompt user to add libs to env path or copy DLLs?
    - Switch pagmo to Git Submodules instead of CMake FetchContent? Possibly better for editing pagmo sources locally. 
    - Pathfinder test
    - Pagmo AVX analysis

# <img src="https://diabetes.zcu.cz/img/icon.png" width="24" height="24" /> SmartCGMS - release repository
This repository is an aggregate repository coupling all parts of the SmartCGMS software framework into a single point. Additionally, every time we release a new version, it becomes available here.

Feel free to fork this repository if you require a customized build setup for your needs.

SmartCGMS software architecture and framework.
Project homepage: [diabetes.zcu.cz](https://diabetes.zcu.cz/smartcgms)

## License

The SmartCGMS software and its components are distributed under the Apache license, version 2. When publishing any derivative work or results obtained using this software, you agree to cite the following paper:

_Tomas Koutny and Martin Ubl_, "SmartCGMS as a Testbed for a Blood-Glucose Level Prediction and/or Control Challenge with (an FDA-Accepted) Diabetic Patient Simulation", Procedia Computer Science, Volume 177, pp. 354-362, 2020

See attached LICENSE file for full licencing information.

|![University of West Bohemia](https://www.zcu.cz/en/assets/logo.svg)|![Department of Computer Science and Engineering](https://www.kiv.zcu.cz/site/documents/verejne/katedra/dokumenty/dcse-logo-barevne.png)|
|--|--|
