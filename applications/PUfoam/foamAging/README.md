@ingroup app_aging

## Aging of polyurethane foam thermal properties
### Physical background
PU foam is blown by CO2 and cyclo-pentane. In time, CO2 is diffusing out of the foam, whereas air is diffusing into the foam. CO2 has lower thermal conductivity than air. Thus, thermal conductivity of the foam is increasing in time.

### Model
Mathematical models in this software simulate the diffusion of gases through the foam and estimate thermal properties of the foam based on its morphology and composition of gases in the cells.

## Installation
LAPACK and BLAS libraries are required. You can install them in Ubuntu using:
```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
```
Compile the detailed models using:
```
./build
```

### (Optional)
Some results can be easily visualized using [VEUSZ](http://home.gna.org/veusz/) program. It can be installed in Ubuntu using:
```
sudo add-apt-repository ppa:jeremysanders/ppa
sudo apt-get update
sudo apt-get install veusz
```

## Aging simulation
First, prepare input file **foamAging.json**. Example input can be found in "example inputs" folder.

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
You have also option to calculate foam conductivity for just one foam at specified conditions using the detailed model. To do this, prepare input file **foamConductivity.json**. Example input can be found in "example inputs" folder.

Run the model using:
```
./foamCond
```
The results can be found in the "results" folder. Look for files named "foamConductivity.out" and "hahtf.out".
