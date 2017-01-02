@ingroup mod_bubbleGrowth

Foam aging
==========

## Scope of this module

This module contains the detailed model for the foam aging simulation. It is a top level application. There is no surrogate model, which needs to be fitted.

## Installation

First make sure you have all dependencies listed in @ref dep_foamAging. A step by step guide is provided in @ref app_aging.

The compilation of the model is done by CMake. You can do
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
