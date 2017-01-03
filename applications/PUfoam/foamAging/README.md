@ingroup app_aging

# Aging of polyurethane foam thermal properties
PU foam is blown by CO<sub>2</sub> and cyclo-pentane. In time, CO<sub>2</sub> is diffusing out of the foam, whereas air is diffusing into the foam. CO<sub>2</sub> has lower thermal conductivity than air. Thus, thermal conductivity of the foam is increasing in time.

Mathematical models in this software simulate the diffusion of gases through the foam and estimate thermal properties of the foam based on its morphology and composition of gases in the cells.

## Installation
### Install LAPACK and BLAS and VEUSZ
You can install them in Ubuntu using:
```
sudo apt-get install liblapack-dev libblas-dev veusz
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

## Running the simulations
1. Prepare input file `foamAging.json` for aging simulation and/or `foamConductivity.json` for the prediction of heat insulation properties. Example input can be found in "examples" folder. The files must be located inside the `inputs` folder. More information about the inputs is provided in `INPUTS.md`.

- Load all surrogate models and their parameters to database:
```
./initModels
```
Note that foamConductivity takes initial points from `init_foamConductivity.json`. This file can be prepared using `prep_init_foamConductivity.py`.
- The heat insulation properties can be predicted by running
```
./workflow_foamConductivity
```
Input morphology and gas composition can be specified manually or taken from results of Bubble growth or 0D or 3D simulation. Note that running this script will delete the parameters of `foamConductivity` from the database. You should run `initModels` again before you run the aging simulation.
- The evolution of heat insulation properties in time can be predicted by running
```
./workflow_foamAging
```
Input morphology and gas composition can be specified manually or taken from results of Bubble growth or 0D or 3D simulation. Note that it is often advantageous to initialize the foamConductivity only for morphology, which would be used in the aging simulation. This can be achieved by running `./prep_init_foamConductivity.py 1`, which will prepare initial points for the initialization based on input file for the aging simulation. However, note that foamConductivity will be re-fitted several times during the simulation, because the gas composition will change and the surrogate model will get out of range. If you don't want any re-fitting during the simulation, then run `./prep_init_foamConductivity.py 2`, which will still use the foam morphology for the initialization of the foamConductivity, but also a whole range of possible gas compositions.

## Results
The results can be found in the `results` folder. They can be visualized using python scripts and VEUSZ.
