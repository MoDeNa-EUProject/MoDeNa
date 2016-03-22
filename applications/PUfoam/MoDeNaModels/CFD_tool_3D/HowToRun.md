# Introduction
This file explains in details how the CFD 3D code should be executed. The test case is a 2-dimensional example of mixing-cup. The recipe has been adopted from Winkler PhD dissertation (Winkler, 2009).  

# Step 1. MoDeNa framework should be installed.
```
cd $HOME
git clone https://github.com/MoDeNa-EUProject/MoDeNa.git
```

### Installing MoDeNa framework:
```
sudo apt-get install automake libltdl-dev libltdl7 mongodb \
    python-rpy2 python-pip python-scipy python-rpy2 python-blessings \
    r-base r-base-dev
sudo apt-get install swig
sudo pip install FireWorks pymongo mongoengine
```
This step should be done within the **R** environment: 
```
cd $HOME
user@machine> R
install.packages('lhs',   repos='http://cran.cnr.Berkeley.edu/', lib=.libPaths()[1], dependencies=TRUE)
install.packages('nlmrt', repos='http://cran.cnr.Berkeley.edu/', lib=.libPaths()[1], dependencies=TRUE)
q()
```
Press `n` to exit **R** without saving the workspace. The next step should be executed from MoDeNa highest directory level. For example, if you have cloned **MoDeNa** in your home directory: 
```
cd $HOME/MoDeNa/src
cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME .
make
make install
```
# Step 2. Install OpenFOAM 2.3.0
In order to install OpenFOAM, you can either follow the steps in [OpenFOAM](http://openfoam.org/archive/2.3.0/download/ubuntu.php) official website or copy/paste the commands below: 
```
cd $HOME
wget http://rheologic.net/sites/default/files/downloads/RheologicRemix-2.3.0-Ubuntu-12.04-LTS.tgz
tar -xvf RheologicRemix-2.3.0-Ubuntu-12.04-LTS.tgz
cat bashrc-OpenFOAM-2.3.0 >> .bashrc 
mkdir -p $FOAM_RUN 
mkdir -p $FOAM_USER_APPBIN
mkdir -p $FOAM_USER_LIBBIN
```
# Step 3. Install MongoDB C++ Driver
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
# Step 4. Compile the source code:
Assuming you have MoDeNa on your HOME directory:
```
cd $HOME/MoDeNa/applications/PUfoam/MoDeNaModels/CFD_tool_3D/src
./makeEigen
wclean
wmake
```
# Step 5. Set the environmental variables
```
export LD_LIBRARY_PATH=/home/mohsen/OpenFOAM/mohsen-2.3.0/platforms/linux64GccDPDebug/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${HOME}/lib/pkgconfig:/usr/local/lib/pkgconfig
export PYTHONPATH=${PYTHONPATH:-}:${HOME}/lib/python2.7/site-packages
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib/python2.7/site-packages:${HOME}/lib/modena:/usr/local/lib
```
# Step 6. Initialize the surrogate models
Since the surrogate models are similar for a 0D or 3D case, you can simply follow the steps on the `README.md`file on: `/MoDeNa/applications/PUfoam/foamingProcess0D`
to initialize the surrogate models.

# Step 7. Run the test case
```
cd $HOME/MoDeNa/applications/PUfoam/MoDeNaModels/CFD_tool_3D/run
./workflow.py
```
### Notes:
-------
* If the following error occurs:
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
* Remember to explicitly initialize the `Kinetics` using:
```
cd $HOME/MoDeNa/publicRepo/MoDeNa/applications/PUfoam/MoDeNaModels/Kinetics
./initModels
```
