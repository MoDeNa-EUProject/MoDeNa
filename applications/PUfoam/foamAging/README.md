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
First, prepare input file **foamAging.in**. Example input can be found in "example inputs" folder. The input file has to have following structure:
- number of time steps, for which the output is produced [-] (integer)
- number of finite volumes for each wall [-] (integer)
- beginning time [s] (float)
- end time [s] (float)
- temperature of aging [K] (float)
- temperature of conductivity measurements [K] (float)
- polymer density [kg/m^3] (float)
- initial pressure [Pa] (float)
- pressure of air, CO2, cyclo-pentane outside the foam [Pa] (float)
- molar fraction of air, CO2, cyclo-pentane in foam [-] (float)
- foam half-thickness [m] (float)
- thickness of foam wall [m] (float)
- cell size [m] (float)
- strut content [-] (float)
- foam density [kg/m^3] (float)
- solubility model (air, CO2, CyP). 0=constant, 1=modena (integer)
- constant solubility. used if solubility model=0. provide any three numbers otherwise (float)
- diffusivity model (air, CO2, CyP). 0=constant, 1=modena (integer)
- constant diffusivity. used if diffusivity model=0. provide any three numbers otherwise (float)

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
You have also option to calculate foam conductivity for just one foam at specified conditions using the detailed model. To do this, prepare input file **foamConductivity.in**. Example input can be found in "example inputs" folder. The input file has to have following structure:
- temperature of boundary 1 - higher [K] (float)
- temperature of boundary 2 - lower [K] (float)
- molar fraction of CO2, air, cyclo-pentane in foam [-] (float)
- emissivity of boundary 1 [-] (float)
- emissivity of boundary 2 [-] (float)
- gas density [kg/m^3] (float)
- polymer density [kg/m^3] (float)
- porosity [-] (float)
- cell size [m] (float)
- morphology input 1=wall thickness, 2=strut content, 3=strut diameter (3 is recommended others can have multiple solutions) [-]
- wall thickness [m] (float)
- strut content [-] (float)
- strut diameter [m] (float)
- foam thickness [m] (float)
- number of finite volumes [-] (integer)
- number of testing rays [-] (integer)
- use wall thickness distribution? - t/f [-] (bool)
- wall thickness standard deviation [-] (float)
- number of gray boxes [-] (integer)
- calcualte effective conductivity numerically? - t/f [-] (bool)
- name of the file with morphology [-] (string)
Run the model using:
```
./foamCond
```
The results can be found in the "results" folder. Look for files named "foamConductivity.out" and "hahtf.out".
