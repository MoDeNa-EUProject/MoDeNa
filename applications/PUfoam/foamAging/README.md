@ingroup app_aging

## Aging of polyurethane foam thermal properties
### Physical background
PU foam is blown by CO2 and cyclo-pentane. In time, CO2 is diffusing out of the foam, whereas air is diffusing into the foam. CO2 has lower thermal conductivity than air. Thus, thermal conductivity of the foam is increasing in time.

### Model
Mathematical models in this software simulate the diffusion of gases through the foam and estimate thermal properties of the foam based on its morphology and composition of gases in the cells.

## Installation
### Install LAPACK and BLAS
You can install them in Ubuntu using:
```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
```
### Install fson library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/japaf/fson.git
cd fson
cmake .
make
make install
```
### Compile the detailed models
```
./build
```

### (Optional)
Some results can be easily visualized using [VEUSZ](http://home.gna.org/veusz/) program. It can be installed in Ubuntu using:
```
sudo apt-get install veusz
```

## Aging simulation
First, prepare input file **foamAging.json**. Example input can be found in "example_inputs" folder. More information about the inputs in INPUTS.md.

Load all surrogate models and their parameters to database:
```
./initModels
```
Main simulation:
```
./workflow
```
The "results" folder contains a file with time dependence of equivalent conductivity "keq_time.out" and files with concentration profiles. The evolution of equivalent conductivity can be viewed using "keq_time.vsz" script and VEUSZ. Concentration profiles can be viewed using "degas_plot.py". Results of the foam conductivity model are in the "launcher" folder.

## Foam conductivity
You have also option to calculate foam conductivity for just one foam at specified conditions using the detailed model. To do this, prepare input file **foamConductivity.json**. Example input can be found in "example_inputs" folder. More information about the inputs in INPUTS.md.

Run the model using:
```
./foamCond
```
The results can be found in the "results" folder. Look for files named "foamConductivity.out" and "hahtf.out".
