@ingroup app_aging

The input files should be located in the `inputs` folder. `foamAging.json` controls aging simulation, whereas `foamConductivity.json` controls foam conductivity simulation. `init_foamConductivity.json` provides initial points for the initialization of foamConductivity.

### Preparing foamAging.json
The following key-pairs should be defined:
- `modelType`: "homogeneous" model first estimates effective foam diffusivity using analytical model, and then simulates diffusion of blowing agents through foam. "heterogeneous" model represents foam as a series of parallel polymer walls separated by gas cells and simulates diffusion of blowing agents directly in this foam. ["homogeneous","heterogeneous"]
- `numerics`:
    - `timeStart`: starting time (s)
    - `timeEnd`: end time (s)
    - `progressTime`: determines spacing between outputs ["logarithmic","linear"]
    - `numberOfOutputs`: how many times should the intermediate results be written, used for linear progress of time
    - `outputsPerOrder`: how many times we save per decade if we have logarithmic progress of time
    - `numberOfOrders`: for how many decades we save in logarithmic progress of time - used to determine total number of outputs, assumes starting time is 0 
    - `wallDiscretization`: number of finite volumes in each wall, used in heterogeneous model
    - `cellDiscretization`: number of finite volumes in each cell, used in heterogeneous model
    - `foamDiscretization`: number of finite volumes in foam, used in homogeneous model
    - `sheetDiscretization`: number of finite volumes in sheet, if sheet is used
- `sourceOfProperty`: when "DirectInput" is used, the property must be given in the input file. Otherwise, it will be loaded from the results of the specified given tool. You can obtain "BubbleGrowth", "Qmom0D", or "Qmom3D" results by running `./workflow_bubbleGrowth`, `./workflow_0D`, or `./workflow_3D` in foamExpansion application. The final state is always saved in `after_foaming.txt` file. "StrutContent" uses StrutContent surrogate model to estimate strut content based on foam density.
    - `foamDensity`: ["DirectInput","BubbleGrowth","Qmom0D","Qmom3D"]
    - `cellSize`: ["DirectInput","BubbleGrowth","Qmom0D","Qmom3D"]
    - `gasComposition`: molar fractions are always normalized ["DirectInput","BubbleGrowth","Qmom0D","Qmom3D"]
    - `strutContent`: ["DirectInput","StrutContent"]
    - `wallThickness`: ["DirectInput"]
- `foamCondition`:
    - `foamHalfThickness`: size of computational domain - half of foam thickness (m)
    - `inProtectiveSheet`: is foam enclosed in a sheet [true,false]
    - `sheetThickness`: thickness of the sheet, if it is used (m)
    - `agingTemperature`: temperature of aging (K)
    - `conductivityTemperature`: temperature of conductivity measurements (K)
    - `initialPressure`: initial pressure in foam (Pa)
    - `initialComposition`: used when "DirectInput" was selected
        - `{gasname}`: molar fraction of gas in initial foam, gases are defined in gasConductivity index set
    - `boundaryPressure`:
        - `{gasname}`: pressure of gas at the outer boundary (Pa)
- `morphology`:
    - `foamDensity`: foam density, if "DirectInput" was selected (kg/m3)
    - `cellSize`: cell size, if "DirectInput" was selected (m)
    - `strutContent`: strut content, if "DirectInput" was selected
    - `wallThickness`: wall thickness, if "DirectInput" was selected (m)
- `physicalProperties`:
    - `polymerDensity`: polymer density (kg/m3)
    - `molarMass`:
        - `{gasname}`: molar mass of gas (kg/mol)
    - `foam`:
        - `solubilityModel`: When "constant" is used, property must be given in the input file. When "modena" is used, surrogate Solubility model is used for given temperature of aging.
            - `{gasname}`: solubility model ["constant","modena"]
        - `solubility`:
            - `{gasname}`: solubility of gas, if "constant" model is used (g/g/bar)
        - `diffusivityModel`: When "constant" is used, diffusivity in polymer must be given in the input file. When "modena" is used, surrogate Diffusivity model is used for given temperature of aging. When "foam" is used, diffusivity in foam (effective diffusivity) must be given in the input file. 
            - `{gasname}`: diffusivity model ["constant","modena","foam"]
        - `diffusivity`:
            - `{gasname}`: diffusivity of gas, if "constant" or "foam" model is used (m2/s)
    - `sheet`:
        - `solubility`:
            - `{gasname}`: solubility of gas in sheet, if "constant" model is used (g/g/bar)
        - `diffusivity`:
            - `{gasname}`: diffusivity of gas in sheet, if "constant" model is used (m2/s)

### Preparing foamConductivity.json
- `upperBoundary`:
    - `emittance`: surface emittance, number 0-1, 0.9 recommended
    - `temperature`: temperature (K)
- `lowerBoundary`:
    - `emittance`: surface emittance, number 0-1, 0.9 recommended
    - `temperature`: temperature (K)
- `gasDensity`: gas density (kg/m3)
- `solidDensity`: polymer density (kg/m3)
- `sourceOfProperty`: when "DirectInput" is used, the property must be given in the input file. Otherwise, it will be loaded from the results of the specified given tool. You can obtain "BubbleGrowth", "Qmom0D", or "Qmom3D" results by running `./workflow_bubbleGrowth`, `./workflow_0D`, or `./workflow_3D` in foamExpansion application. The final state is always saved in `after_foaming.txt` file. "StrutContent" uses StrutContent surrogate model to estimate strut content based on foam density.
    - `porosity`: ["DirectInput","BubbleGrowth","Qmom0D","Qmom3D"]
    - `cellSize`: ["DirectInput","BubbleGrowth","Qmom0D","Qmom3D"]
    - `gasComposition`: molar fractions are always normalized ["DirectInput","BubbleGrowth","Qmom0D","Qmom3D"]
    - `strutContent`: ["DirectInput","StrutContent"]
    - `wallThickness`: ["DirectInput"]
- `gasComposition`:
    - `{gasname}`: molar fraction of gas, gases are defined in gasConductivity index set
- `porosity`: foam porosity, if "DirectInput" is used
- `cellSize`: cell size, if "DirectInput" is used (m)
- `morphologyInput`: when we know foam porosity and cell size, foam morphology is fully specified, if we additionally know either wall thickness, strut content, or strut size. Specify, which property you know, the others will be used as initial guesses to calculate them precisely. "strutSize" is recommended. "strutSize2" uses slightly different model to calculate strut size and wall thickness. ["wallThickness","strutContent","strutSize","strutContent2"]
- `strutContent`: strut content, if "DirectInput" is used
- `wallThickness`: wall thickness, if "DirectInput" is used (m)
- `strutSize`: strut size (m)
- `foamThickness`: foam thickness (m)
- `spectra`:
    - `polymer_n`: file with spectrum of real part of refractive index of polymer
    - `polymer_k`: file with spectrum of imaginary part of refractive index of polymer
    - `gas_k`: file with spectrum of imaginary part of refractive index of gas
- `spatialDiscretization`: number of finite volumes for conduction-radiation simulation
- `useWallThicknessDistribution`: assume wall thickness distribution for calculation of radiative properties [true,false]
- `wallThicknessStandardDeviation`: if wall thickness distribution is used
- `numberOfGrayBoxes`: number of gray boxes for discretization of foam spectra, 10 is recommended
- `numericalEffectiveConductivity`: use numerical instead of analytical model to calculate effective conductivity, false is recommended [true,false]
- `structureName`: filename with foam morphology saved in VTK voxel format, if numerical model for effective conductivity is used
- `testMode`: disable estimation of radiative properties, significantly decreases computational time, false is recommended [true,false]

### Preparing init_foamConductivity.json
- `T`: list of temperatures (K)
- `dcell`: list of cell sizes (m)
- `eps`: list of porosities
- `fstrut`: list of strut contents
- `x[{gasname}]`: list of gas molar fractions, gases are defined in gasConductivity index set

### Preparing cell_gas.json
- `temperature`: 
    - `min`: minimum of temperature interval (K),
    - `max`: maximum of temperature interval (K),
    - `points`: number of points in temperature interval
- `initial_weight_fraction`: 
    `{gasname}`: initial weight fraction in foam recipe, blowing agents defined in gasConductivity index set and H2O are supported
- `molar_mass`: 
    `{gasname}`: molar mass for H2O and all used gases (kg/mol)
- `polymer_density`: polymer density (kg/m3),
- `foam_density`: foam density (kg/m3),
- `cell_size`: cell size (m)