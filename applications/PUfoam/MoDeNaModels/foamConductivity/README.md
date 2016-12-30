@ingroup mod_foamConductivity

Foam aging
==========

## Scope of this module

This module contains the detailed model for the foam conductivity simulation. It can be run as a top level application. However, it also contains a surrogate model, which can be embedded into other models.

## Installation

First make sure you have all dependencies listed in @ref dep_foamConductivity. A step by step guide is provided in @ref app_aging.

The compilation is done by CMake. You can do
```
cd src
cmake .
make
```
However, you can also use the build script in @ref app_aging.

## Input / Output

The details on input and output are provided in @ref app_aging.

## Run

Although you can start the executable directly, it is strongly advised that you use the workflow script in @ref app_aging.
