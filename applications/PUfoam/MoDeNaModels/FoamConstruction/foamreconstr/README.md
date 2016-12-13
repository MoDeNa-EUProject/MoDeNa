# Program for foam reconstruction
- based on Voronoi tessellation
- closed-cell foams
- open-cell foams
- struts
- voxel based morphology

## Installation
You must first install voro++ library http://math.lbl.gov/voro++/. You can try
to execute following script, which should download, compile and install it:
```
sudo ./install_voro++.sh
```
By default it is installed in `/usr/local`. You must change `CMakeLists.txt`
in `src/` if you installed it somewhere else.
The program depends on `GSL` library. Make sure it is installed by running:
```
sudo apt-get install gsl-bin libgsl0-dev
```
The program is compiled by:
```
cmake .
make
```

## Set parameters
Edit `foamreconstr.in`:
1. Create nodes? `bool`
1. Create edges? `bool`
1. Open-cell foam? `bool`
1. Save for DX? `bool`
1. Save for Paraview? `bool`
1. How large the nodes should be. `double`
1. How large the edges should be. `double`
1. Desired porosity of only struts. `double`
1. Perturbation of seed positions. `double` in range 0..1
1. General seed positions. 1=cubic grid, 2=hexagonal lattice, ABAB...,3=hexagonal lattice, ABCABC...
1. Size of the domain in `x` in voxels. `int`
1. Size of the domain in `y` in voxels. `int`
1. Size of the domain in `z` in voxels. `int`
1. Size of the cells in `x` in voxels. `int`
1. Size of the cells in `y` in voxels. `int`
1. Size of the cells in `z` in voxels. `int`
1. Keep copy of tessellation - gnuplot diagram? `bool`
1. Keep copy of tessellation - alternative gnuplot diagram? `bool`
1. Import morphology from file? In this case domain size is set from input file. `bool`
1. Display detailed progress report to standard output. `bool`
1. Name of output file with morphology without extension. `string`
1. Name of input VTK file with morphology. `string`
1. Name of input/output file with tessellation - gnuplot diagram. `string`
1. Name of output file with tessellation - alternative gnuplot diagram. `string`
1. Name of output file with morphology descriptors. `string`
1. Name of output file with parameters. `string`

## Run
```
./foamreconstr
```

## Postprocess
`*.vtk` file can be viewed in Paraview. File with morphology descriptors
contains porosity and strut content.
