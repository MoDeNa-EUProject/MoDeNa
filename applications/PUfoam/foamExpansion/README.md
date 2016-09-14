@ingroup app_foaming

# Foaming simulation
Collection of software tools for the simulation of polyurethane (PU) foaming
process. Based on a recipe the model predicts the evolution of chemical
kinetics, temperature, foam density, bubble size distribution, etc.
This directory includes the foaming process in zero dimensional and three
dimensional spaces. Additionally, a wall drainage simulation can be performed, which estimates the size and shape of strut and wall thickness profile. In order to run a test case, firstly the source codes should
be compiled. The script `build` compiles the necessary models. The surrogate
models will be loaded into the database by executing the two scripts provided:
`initModels` and `initBubbleGrowth`. Finally, the scripts `workflow*` run the
simulations and detailed model for bubble growth.

# Dependencies and Installation
Several models need to be compiled. C++ and Fortran compilers are required.
First the MoDeNa framework should be compiled.
The [following steps](https://github.com/MoDeNa-EUProject/MoDeNa) describe how
to install MoDeNa framework.

### 1. Install OpenFOAM
The CFD solver is compiled with OpenFOAM 4.0 on Ubuntu 16.04. In order to install OpenFOAM, you can either follow the steps in [OpenFOAM](http://openfoam.org/download/4-0-ubuntu/) website or
copy/paste the commands below which install OpenFOAM and ParaView:

```
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo sh -c "wget -O - http://dl.openfoam.org/gpg.key | apt-key add -"
sudo apt-get update
sudo apt-get -y install openfoam4
mkdir -p $FOAM_RUN
mkdir -p $FOAM_USER_APPBIN
mkdir -p $FOAM_USER_LIBBIN
```
*Note 1*: This [video](https://www.youtube.com/watch?v=x6BWTEcuh38) demonstrates the installation process.

*Note 2*: Consult with [official webpage](http://openfoam.org/download/4-0-ubuntu/) for installation problems.

### 2. Install MongoDB C++ Driver
```
cd $HOME
git clone -b r1.3 https://github.com/mongodb/mongo-c-driver
cd mongo-c-driver
./autogen.sh
make
sudo make install
cd $HOME
git clone -b master https://github.com/mongodb/mongo-cxx-driver
cd mongo-cxx-driver/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ..
sudo make
sudo make install
```
### 3. Install boost, lapack and blas library
```
sudo apt-get install libboost-dev liblapack-dev libblas-dev
```
### 4. Install PETSc globally:
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
### 5. Install rapidjson library
```
cd where-you-want-source-files
git clone https://github.com/miloyip/rapidjson.git
cd rapidjson
cmake .
make
sudo make install
```
### 6. Install fson library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/japaf/fson.git
cd fson
cmake .
make
make install
```
### 7. Install bspline library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/japaf/bspline-fortran.git
cd bspline-fortran
cmake .
make
make install
```
### 8. Install sundials library
```
cd where-you-want-source-files
git clone https://github.com/luca-heltai/sundials.git
cd sundials
mkdir build
cd build
cmake -DFCMIX_ENABLE=ON -DLAPACK_ENABLE=ON ..
make
sudo make install
```
### 9. Set the environmental variables for MoDeNa
```
user=$(whoami)
export LD_LIBRARY_PATH=/home/${user}/OpenFOAM/${user}-2.3.0/platforms/linux64GccDPDebug/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${HOME}/lib/pkgconfig:/usr/local/lib/pkgconfig
export PYTHONPATH=${PYTHONPATH:-}:${HOME}/lib/python2.7/site-packages
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib/python2.7/site-packages:${HOME}/lib/modena:/usr/local/lib

```
### 10. Compile the models
```
./build
```
Please note that the output of build should not contain any error messages.
Otherwise, the tool cannot work properly.
### 11. Intall VEUSZ plotting program (Optional)
Some results can be easily visualized using [VEUSZ](http://home.gna.org/veusz/)
program. It can be installed in Ubuntu using:
```
sudo apt-get install veusz
```

# Run
In order to execute the program, you need to first provide the input variables
for all the modelling tools. In the `examples` directory several cases have
been provided. One can copy the `inputs` directory into the `foamExpansion` and
modify it for the targeted recipe. The `unifiedInput.json` provides the input
data for the zero dimensional simulation, `wallDrainage_inputs.json` provides the inputs for wall drainage simulation and the rest of input directory
creates the test case for OpenFOAM simulation. The details of the input
variables have been elaborated in `INPUTS.md`. After preparing the inputs the following
steps should be executed:

1. Initialize surrogate nano-scale models:
```
./initModels1
./initModels2
```
2.  Execute meso-scopic bubble growth simulation:
```
./workflow_bubbleGrowth
```
3. Initialize surrogate bubble growth model:
```
./initBubbleGrowth
```
4. Run macroscopic simulation (0D or 3D):
```
./workflow_0D
```
or
```
./workflow_3D
```
    *Note:*

    - In case of any changes in the input files steps 2-4 should repeated.

    - Initial moments of the bubble size distribution can be calculated by running
    `initMoments`. This script uses the mean, variance and initial number
    density of bubbles from `unifiedInput.json`.

    - If the following error occurs:
    ```
         --> FOAM FATAL ERROR:
         Wrong number of arguments, expected 0 found 1

         FOAM exiting

         cannot find system Renviron
         Fatal error: unable to open the base package
    ```
    The workaround is to set the R environment variable as below:
    ```
    R_HOME=/usr/lib/R (where you have installed R)
    export R_HOME=/usr/lib/R
    ```

- The wall drainage simulation can be executed anytime after running the `./workflow_bubbleGrowth` by:
```
./workflow_wallDrainage
```

# Results
Results of the last simulation are stored in the corresponding sub-directory of
the `results` directory. As mentioned before, we use `VEUSZ` to visualize the
results of simulations. For example, to look at selected results for a 0D
simulation, you can open `plotQmom0D.vsz` using VEUSZ. Further, the results of
the 3D simulation (stored in the launcher directory) can be displayed using a
third party software such as [paraview.](http://www.paraview.org/) Additional results of wall drainage simulation can be displayed by running `./plotWallDrainage.py`.
