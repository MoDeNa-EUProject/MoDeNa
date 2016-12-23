@ingroup mod_bubbleGrowth

Bubble growth
=============

## Scope of this module

This module contains the detailed model for the bubble growth simulation. A little bit unusually, the surrogate model is created by the initBubbleGrowth script in @ref app_foaming. The description of the model is provided in \cite Ferkl2016.

## Installation

First make sure you have all dependencies listed in @ref dep_bubbleGrowth. A step by step guide is provided in @ref app_foaming.

The compilation is done by CMake. You can do
```
cd src
cmake .
make
```
However, you can also use the build script in @ref app_foaming.

## Input / Output

The details on input and output are provided in @ref app_foaming.

## Run

Although you can start the executable directly, it is strongly advised that you use the workflow script in @ref app_foaming.
