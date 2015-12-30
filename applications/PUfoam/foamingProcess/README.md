@ingroup app_foaming

## Foaming simulation
Collection of software tools for the simulation of polyurethane (PU) foaming
process. Based on a recipe the model predicts the evolution of chemical
kinetics, temperature, foam density, bubble size distribution, etc.

## Installation
Several models need to be compiled. C++ and Fortran compilers are required.
Execute following command:
```
./build
```

### (Optional)
Some results can be easily visualized using [VEUSZ](http://home.gna.org/veusz/)
program. It can be installed in Ubuntu using:
```
sudo add-apt-repository ppa:jeremysanders/ppa
sudo apt-get update
sudo apt-get install veusz
```

## Run
First, you need to specify recipe. Copy/create/modify **unifiedInput.json**. An
example recipies can be found in "example inputs" folder. Create inputs to
models:
```
./createInputs.py
```
Initialize surrogate nano-scale models:
```
./initModels
```
Execute meso-scopic bubble growth simulation:
```
./workflow1
```
Initialize surrogate bubble growth model:
```
./initBubbleGrowth
```
Run macroscopic simulation:
```
./workflow2
```

The above sequence can be performed by executing
```
./run
```

## Results
Results are stored in the "results" folder. To look at selected macroscopic
results you can use "plotQmom.vsz" and VEUSZ.
