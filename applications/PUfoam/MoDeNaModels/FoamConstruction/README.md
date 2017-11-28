@ingroup mod_foamConstruction

Foam Construction
=================

## Scope of this module
This module contains an utility, which can create a spatially three-dimensional image of foam morphology with desired foam density, cell size distribution and strut content. This module does not contain any MoDeNa model.

## Installation
The code depends on several third-party applications. They are listed in @ref dep_foamConstruction. 

To install all of these on Ubuntu, do:
```
sudo apt-get install libmatheval-dev gmsh gsl-bin libgsl0-dev python-vtk \
    lib3ds-1-3 libjpeg62 freeglut3 libnlopt-dev libboost-dev \ libboost-date-time-dev libboost-thread-dev zlib1g-dev libpng12-dev \
    libjpeg8-dev libtiff5-dev libopenexr-dev
```
Install python packages using pip
```
sudo -H pip install spack vapory
```
Install packing-generation from github and copy the executable to PATH
```
cd where-you-want-source-files
git clone https://github.com/VasiliBaranov/packing-generation.git
cd packing-generation/_Release
make
cp PackingGeneration.exe ~/bin/
```
Then download and install `neper` version 3 from http://neper.sourceforge.net/downloads.html.
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
- `packing.py` - sphere packing
- `tessellation.py` - tessellation, creation of RVE
- `periodicBox.py` - creation of the box with periodic boundary conditions
- `vtkconv.py` - conversion from binary vtk to ascii vtk
- `geo_tools.py` - manipulation of GMSH geometry files, creation of unstructured mesh
- `run.py` - main executable script
- `simulation.py` - experimental script for FEM simulations using FEniCS (conductivity or diffusivity)

All files must be in one directory.

## Inputs
The code is controlled by the `input.json` file, which must be located in the
root of `FoamConstruction` folder. Default input file can be found
in `example_inputs` directory. Following inputs can be adjusted:
- `filename`: base name of created files 
- `packing`: create sphere packing [true, false], 
- `packing_options`: 
    - `shape`: shape of log-normal distribution, 
    - `domain_size`: domain size (1 is recommended), 
    - `number_of_cells`: number of cells, 
    - `scale`: scale of log-normal distribution, 
    - `algorithm`: type of sphere packing algorithm [`simple`, `-ls`, `-fba`, `-lsgd`, `-lsebc`, `-ojt`, `-kjt`], `simple` algorithm is included, others are enabled by packing-generation (see https://github.com/VasiliBaranov/packing-generation), `-fba` is recommended
- `tessellation`: create tessellated foam [true, false], 
- `tessellation_options`: 
    - `visualize_tessellation`: visualize tessellation [true, false], false is recommended
- `structured_grid`: create structured (voxel) mesh [true, false],
- `structured_grid_options`: 
    - `render_box`: visualize foam [true, false], false is recommended 
    - `strut_content`: strut content, 
    - `porosity`: foam porosity, 
    - `strut_size_guess`: strut size in voxels, guess usually 4-8
    - `binarize_box`: run part of the script, which creates voxel mesh, [true, false], true is recommended
    - `move_to_periodic_box`: run part of the script, which moves foam to periodic box, [true, false], true is recommended
- `unstructured_grid`: create unstructured (tetrahedral) mesh [true, false], 
- `unstructured_grid_options`: 
    - `create_geometry`: run part of the script, which creates foam, [true, false], true is recommended, 
    - `convert_mesh`: run part of the script, which converts mesh to .xml, [true, false], true is recommended, 
    - `wall_thickness`: wall thickness parameter, 0.02 is good guess
    - `mesh_domain`: run part of the script, which creates mesh, [true, false], true is recommended

## Execution
Prepare `input.json`, then:
```
./run.py
```
Optimizing porosity and strut content of voxelized foam is relatively time
consuming. You can switch to `Bounded` method if you approximately know the
size of the box in voxels (usually from experience with the program). In that
case you need to edit the `run.py` script.

## Outputs
Several output files are created. Structured mesh is in `{filename}_str.vtk`, unstructured mesh is in `{filename}_uns.vtk`.

Generally, `.geo`, and `.msh` files can be viewed with `gmsh`. `.stl`, `.ply` and `.vtk`
files can be viewed with `paraview`.
