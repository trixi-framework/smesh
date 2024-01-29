# smesh

[![Build Status](https://github.com/trixi-framework/smesh/workflows/CI/badge.svg)](https://github.com/trixi-framework/smesh/actions?query=workflow%3ACI)
[![Coveralls](https://coveralls.io/repos/github/trixi-framework/smesh/badge.svg)](https://coveralls.io/github/trixi-framework/smesh)
[![Codecov](https://codecov.io/gh/trixi-framework/smesh/branch/main/graph/badge.svg)](https://codecov.io/gh/trixi-framework/smesh)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/license/mit/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10579422.svg)](https://doi.org/10.5281/zenodo.10579422)

A simple Fortran package for generating and handling unstructured triangular and polygonal
meshes.


## Getting started
### Prerequisites
* CMake v3.5.1
* a somewhat recent Fortran compiler
  * tested with gfortran v11 Linux, macOS
  * tested with gfortran v13 on Linux, macOS, Windows

### Installation
To use smesh, you need to compile it first. We test the compilation regularly using our CI
setup with gfortran on Linux, macOS, and Windows (the latter via MSYS2).

To build and install, perform the following steps:
* Get the sources (e.g., by cloning this repository)
* Create a `build` directory for intermediate build artifacts
* Configure with CMake
* Build the library and executable products
* Install everything

On most systems, the following commands should achieve to build and install smesh into a
local directory:
```shell
git clone git@github.com:trixi-framework/smesh.git
# Alternatively use this if you do not have your system set up for git via ssh for GitHub:
# git clone https://github.com/trixi-framework/smesh.git
mkdir smesh/build && cd smesh/build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
cmake --build .
cmake --install .
cd ..
```
This will install everything into the `smesh/install` directory.

### Usage
To use smesh, you can call the `smesh_run` executable with the path to a smesh-compatible
configuration file as the first command line argument, e.g.,
```shell
cd smesh
install/bin/smesh_run smesh_example.cfg
```
This will give you an output similar to
```
Computing Delaunay triangulation.
Triangulation elements:        775
Total flipped edges:          1248
Average search time:          6.68
Flips/triangle:               1.61
Flips/node:                   3.09
```
and some additional output files `*.dat` in the current directory.


## Referencing
If you use smesh in your own research, please cite this repository as follows:
```bibtex
@misc{chiocchetti2024smesh,
  title={smesh: {A} simple {F}ortran package for generating and handling unstructured triangular and polygonal meshes},
  author={Chiocchetti, Simone},
  year={2024},
  howpublished={\url{https://github.com/trixi-framework/smesh}},
  doi={10.5281/zenodo.10579422}
}
```


## Authors
Smesh was initiated by
[Simone Chiocchetti](https://www.mi.uni-koeln.de/NumSim/dr-simone-chiocchetti/)
(University of Cologne, Germany), who is also its principal maintainer.


## License and contributing
Smesh is available under the MIT license (see [LICENSE.md](LICENSE.md)).
Contributions by the community are very welcome!
