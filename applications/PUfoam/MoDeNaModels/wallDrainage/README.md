@ingroup mod_wallDrainage

Wall drainage
=============

## Scope of this module

This module contains the detailed model for the wall drainage simulation. It is run as a top level application. It can use results previously calculated by Bubble growth or CFD_0D or CFD_3D models.

## Installation

First make sure you have all dependencies listed in @ref dep_wallDrainage. A step by step guide is provided in @ref app_foaming, but necessary steps are also provided below.

### 1. Install fson library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/josephalevin/fson.git
cd fson
make
make install
```
### 2. Install bspline library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/japaf/bspline-fortran.git
cd bspline-fortran
cmake .
make
make install
```

The compilation of the model is done by CMake. You can do
```
cmake .
make
```

## Input / Output
The details on input and output are provided in @ref app_foaming. You can visualize results using plot.py.

## Run
Although you can start the executable directly, it is strongly advised that you use the workflow script in @ref app_foaming. The direct execution is performed by
```
./walldrainage
```
