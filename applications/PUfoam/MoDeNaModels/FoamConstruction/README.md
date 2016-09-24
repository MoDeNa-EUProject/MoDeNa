Foam Construction
============================
## Installation
The code depends on several third-party applications:
- `neper` and `voro++` for tessellation
- `POV-Ray` for visualization of tessellation
- `gmsh`, `vtk` and `meshconv` for mesh manipulation
- `binvox` for voxelization
- `foamreconstr` for creation of voxelized foams with struts

To install all of these on Ubuntu, do:
```
sudo apt-get install libmatheval-dev gmsh gsl-bin libgsl0-dev python-vtk \
    lib3ds-1-3 libjpeg62 freeglut3 libnlopt-dev libboost-dev \ libboost-date-time-dev libboost-thread-dev zlib1g-dev libpng12-dev \
    libjpeg8-dev libtiff5-dev libopenexr-dev
```
Then download and install `neper` from http://neper.sourceforge.net/downloads.html.
You will need to unpack `neper`, go to its `src` folder and then:
```
mkdir build
cd build
cmake ..
make
sudo make install
```
Install POV-Ray from sources:
```
cd where-you-want-source-files
git clone https://github.com/POV-Ray/povray.git
cd povray/unix
./prebuild.sh
cd ..
./configure COMPILED_BY="your name <email@address>"
make
sudo make install
```
Download `meshconv` from http://www.patrickmin.com/meshconv/, download
`binvox` from http://www.patrickmin.com/binvox/ and copy `meshconv` and `binvox` to `$PATH`. You can do this manually or in terminal using
```
wget http://www.patrickmin.com/meshconv/linux64/meshconv
chmod +x meshconv
sudo mv meshconv /usr/local/bin/meshconv
wget http://www.patrickmin.com/binvox/linux64/binvox
chmod +x binvox
sudo mv binvox /usr/local/bin/binvox
```
To install `voro++`:
```
cd foamreconstr
sudo ./install_voro++.sh
```
If that doesn't work, you will need to download it from
http://math.lbl.gov/voro++/download/ and install it manually.

To compile `foamreconstr` go to `foamreconstr/` folder and:
```
cmake .
make
```

## Files
The folder `FoamConstruction` contains following files:
- `FoamGeometryConstruction_Periodic.py` - tessellation, creation of RVE
- `periodicBox.py` - creation of the box with periodic boundary conditions
- `vtkconv.py` - conversion from binary vtk to ascii vtk
- `run` - main executable script

All files must be in one directory.

## Inputs
The code is controlled by the `input.json` file, which must be located in the
root of `FoamConstruction` folder. Default input file can be found
in `example_inputs` directory. Following inputs can be adjusted:
- `MU` - mean of cell size distribution
- `SIGMA` - standard deviation of cell size distribution
- `NumOfCells` - number of cells in RVE (representative volume element)
- `porosity` - desired porosity of voxelized foam
- `strutContent` - desired strut content of voxelized foam
- `filename` - name of the output file with RVE
- `deleteFiles` - delete some redundant output files after execution
- `packing` - call packing algorithm, which creates seeds and radii for
tessellation
- `alternativePackingAlgorithm` - uses alternative packing algorithm from "spherepack"
- `tesselation` - call tessellation program, which creates foam with desired
cell size distribution based on results of packing
- `visualizeTesselation` - visualizes tessellation using the POV-Ray
- `geometry` - limits the foam to RVE
- `statistics` - compute and save cell volumes, face surfaces, etc.
- `hypermesh` - create input for `Hypermesh`
- `moveToPeriodicBox` - move RVE to box with periodic boundary conditions
- `renderBox` - show the box with periodic boundary conditions
- `binarizeBox` - voxelize the box with periodic boundary conditions

## Execution
Prepare `input.json`, then:
```
./run.py
```
Optimizing porosity and strut content of voxelized foam is relatively time
consuming. You can switch to `Bounded` method if you approximately know the
size of the box in voxels (usually form experience with the program). In that
case you need to edit the `run` script.

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
- `PeriodicRVEBox.ply` - surface mesh of the box with periodic boundary
conditions
- `PeriodicRVEBox.vtk` - voxelized version of the box with periodic
boundary conditions - main output for foam with no struts
- `PeriodicRVEBox-ascii.vtk` - voxelized version of the box with periodic
boundary conditions - ascii version
- `PeriodicRVEBoxStruts.vtk` - voxelized version of the box with periodic
boundary conditions - main output for foam with struts

Generally, `.geo` files can be viewed with `gmsh`, `.stl`, `.ply` and `.vtk`
files can be viewed with `paraview`.
