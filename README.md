# Commands

- Git submodules

`git submodule update --remote --init`

- CMake configure step

`cmake --preset msvc-ninja-configure`

- CMake build step

`cmake --build --preset msvc-ninja-build`

- CMake configure + build

`cmake --workflow --preset msvc-ninja-all`

- Set your MSVC and Qt path in CMakePresets.json before attempting to configure or build. 
- To compile outside of IDE use x64 Native Tools Command Prompt for VS

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
