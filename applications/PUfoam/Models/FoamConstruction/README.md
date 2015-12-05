Foam Construction
============================
## Installation
The code depends on several third-party applications:
- `wine` to have ability of running windows executable files (i.e. file.exe).
- `neper` for tessellation
- `gmsh`, `vtk`, `meshconv` and `binvox` for mesh manipulation

To install all of these on Ubuntu, do:
```
sudo add-apt-repository ppa:ubuntu-wine/ppa
sudo apt-get update
sudo apt-get install wine1.7 libmatheval-dev gmsh gsl-bin libgsl0-dev \
    python-vtk
```
Then download and install `neper` from http://neper.sourceforge.net/downloads.html
(follow its README for instructions),
download `meshconv` from http://www.cs.princeton.edu/~min/meshconv/, download
`binvox` from http://www.cs.princeton.edu/~min/binvox/ and copy `meshconv` and
 `binvox` to `$PATH`

## Files
The folder `FoamConstruction` contains following files:
- `SpherePackFB.exe` - packing algorithm
- `FoamGeometryConstruction_Periodic.py` - tessellation, creation of RVE
- `periodicBox.py` - creation of the box with periodic boundary conditions
- `vtkconv.py` - conversion from binary vtk to ascii vtk
- `run` - main executable script

All files must be in one directory.

## Inputs
The code is controlled by the `input.json` file. Default input can be found
in `example_inputs` directory. Following inputs can be adjusted:
- `MU` - mean of cell size distribution
- `SIGMA` - standard deviation of cell size distribution
- `NumOfCells` - number of cells in RVE (representative volume element)
- `porosity` - desired porosity of voxelized foam
- `filename` - name of the output file with RVE
- `deleteFiles` - delete some redundant output files after execution
- `packing` - call packing algorithm, which creates seeds and radii for
tessellation
- `tesselation` - call tessellation program, which creates foam with desired
cell size distribution based on results of packing
- `geometry` - limits the foam to RVE
- `statistics` - compute and save cell volumes, face surfaces, etc.
- `hypermesh` - create input for `Hypermesh`
- `moveToPeriodicBox` - move RVE to box with periodic boundary conditions
- `renderBox` - show the box with periodic boundary conditions
- `binarizeBox` - voxelize the box with periodic boundary conditions

## Execution
Prepare `input.json`, then:
```
./run
```

## Outputs
Several output files are created:
- `PeriodicRVE.geo` which is the periodic geometry of foam in .GEO
format.
- `CellVolumes_PeriodicRVE.txt` containing the volumes of all the
cells in RVE.
- `FaceAreas_PeriodicRVE.txt` containing the areas of all the faces
in RVE.
- `EdgeLengths_PeriodicRVE.txt` containing the lengths of all the
edges in RVE.
- `HyperMeshinput.cmf` which is a input file for `Hypermesh` to
recreate the geometry.
- `PeriodicRVEBox.stl` - surface mesh of the box with periodic boundary
conditions
- `PeriodicRVEBox-ascii.vtk` - voxelized verion of the box with periodic
boundary conditions

Generally, `.geo` files can be viewed with `gmsh`, `.stl`, `.ply` and `.vtk`
files can be viewed with `paraview`.
