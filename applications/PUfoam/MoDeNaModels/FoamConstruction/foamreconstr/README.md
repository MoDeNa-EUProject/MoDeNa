@ingroup src_mod_foamConstruction

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
- Create nodes? `bool`
- Create edges? `bool`
- Open-cell foam? `bool`
- Save for DX? `bool`
- Save for Paraview? `bool`
- How large the nodes should be. `double`
- How large the edges should be. `double`
- Desired porosity of only struts. `double`
- Perturbation of seed positions. `double` in range 0..1
- General seed positions. 1=cubic grid, 2=hexagonal lattice, ABAB...,3=hexagonal lattice, ABCABC...
- Size of the domain in `x` in voxels. `int`
- Size of the domain in `y` in voxels. `int`
- Size of the domain in `z` in voxels. `int`
- Size of the cells in `x` in voxels. `int`
- Size of the cells in `y` in voxels. `int`
- Size of the cells in `z` in voxels. `int`
- Keep copy of tessellation - gnuplot diagram? `bool`
- Keep copy of tessellation - alternative gnuplot diagram? `bool`
- Import morphology from file? In this case domain size is set from input file. `bool`
- Display detailed progress report to standard output. `bool`
- Name of output file with morphology without extension. `string`
- Name of input VTK file with morphology. `string`
- Name of input/output file with tessellation - gnuplot diagram. `string`
- Name of output file with tessellation - alternative gnuplot diagram. `string`
- Name of output file with morphology descriptors. `string`
- Name of output file with parameters. `string`

## Run
```
./foamreconstr
```

## Postprocess
`*.vtk` file can be viewed in Paraview. File with morphology descriptors
contains porosity and strut content.
