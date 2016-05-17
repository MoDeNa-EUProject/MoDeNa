@ingroup app_foaming

## Foaming simulation
Collection of software tools for the simulation of polyurethane (PU) foaming
process. Based on a recipe the model predicts the evolution of chemical
kinetics, temperature, foam density, bubble size distribution, etc.

## Installation
Several models need to be compiled. C++ and Fortran compilers are required.

Install boost library
```
sudo apt-get install libboost-dev
```
Install PETSc globally:
```
sudo mkdir /opt/petsc
sudo chown user:group /opt/petsc
cd /opt/petsc
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.5.tar.gz
tar -xzf petsc-3.4.5.tar.gz
cd petsc-3.4.5
./configure --with-cc=gcc --with-fc=gfortran --download-f-blas-lapack --download-mpich --download-scalapak=yes
make
```
Export variables in `.bashrc`
```
export PETSC_DIR=/opt/petsc/petsc-3.4.5
export PETSC_ARCH=arch-linux2-c-debug
```
Install rapidjson library
```
cd where-you-want-source-files
git clone https://github.com/miloyip/rapidjson.git
cd rapidjson
cmake .
make
sudo make install
```
Compile the models using:
```
./build
```

### (Optional)
Some results can be easily visualized using [VEUSZ](http://home.gna.org/veusz/)
program. It can be installed in Ubuntu using:
```
sudo apt-get install veusz
```

## Run
First, you need to specify recipe. Copy/create/modify **unifiedInput.json**. An
example recipies can be found in "examples" folder.
Initialize surrogate nano-scale models:
```
./initModels
```
Execute meso-scopic bubble growth simulation:
```
./workflow_bubbleGrowth
```
Initialize surrogate bubble growth model:
```
./initBubbleGrowth
```
Run macroscopic simulation:
```
./workflow_0D
```

## Results
Results are stored in the "results" folder. To look at selected macroscopic
results you can use "plotQmom.vsz" and VEUSZ.
