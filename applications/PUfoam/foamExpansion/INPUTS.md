@ingroup app_foaming

The recipe for the simulation is specified in the `inputs` directory. Several test cases are prepared in the `examples` directory. The `unifiedInput.json` file controls the Bubble growth model and the 0D simulation, , `wallDrainage_inputs.json` controls the Wall drainage simulation and the rest of the files and directories control the OpenFOAM simulation.

### Preparing unifiedInput.json
The following key-pairs should be defined:
- `physicalBlowingAgent`: name of blowing agent [`n-pentane`, `R11`]
- `initialConditions`:
    - `concentrations`: initial concentrations of reactants and blowing agents (mol/m^3) (always water, PBL, CO2 and for `Baser` kinetics: isocyanate, polyol; for `RF-1` kinetics: catalyst, polyol1, polyol2, amine, isocyanate1, isocyanate2, isocyanate3)
    - `bubbleRadius`: initial bubble radius (m),
    - `bubbleRadiusDeviation`: standard deviation of initial BSD,
    - `numberBubbleDensity`: initial number density of bubbles (1/m^3),
    - `temperature`: initial temperature (K)
- `kinetics`:
    - `kineticModel`: name of kinetic model [`Baser`, `RF-1`],
    - `useDilution`: use correction for dilution in Baser model (t/f),
    - `gelPoint`: conversion at the gel point
    - `gellingReaction`: parameters of gelling reaction
        - `frequentialFactor`: pre-exponential factor for Baser model (m^3/mol/s),
        - `activationEnergy`: activation energy for Baser model (J/mol),
        - `reactionEnthalpy`: reaction enthalpy used in all models (J/mol)
    - `blowingReaction`:
        - `frequentialFactor`: pre-exponential factor for Baser model (1/s),
        - `activationEnergy`: activation energy for Baser model (J/mol),
        - `reactionEnthalpy`: reaction enthalpy used in all models (J/mol)
- `physicalProperties`:
    - `pressure`: (Pa)
    - `blowingAgents`:
        - `PBL`:
            - `molarMass`: (kg/mol)
            - `heatCapacityInLiquidPhase`: (J/kg/K)
            - `heatCapacityInGaseousPhase`: (J/kg/K)
            - `evaporationHeat`:  (J/kg)
            - `density`: (kg/m^3)
            - `diffusivityModel`: [`constant`,`nanotools`],
            - `diffusivity`: diffusivity when `constant` diffusivityModel is used (m^2/s)
            - `solubilityModel`: [`constant`, `pcsaft`, `Gupta`, `Winkler`, `Baser`] Gupta and Winkler are models for n-pentane, Baser is model for R11
        - `CO2`:
            - `molarMass`: (kg/mol)
            - `heatCapacityInLiquidPhase`: (J/kg/K)
            - `heatCapacityInGaseousPhase`: (J/kg/K)
            - `evaporationHeat`:  (J/kg)
            - `density`: (kg/m^3)
            - `diffusivityModel`: [`constant`,`nanotools`],
            - `diffusivity`: diffusivity when `constant` diffusivityModel is used (m^2/s)
            - `solubilityModel`: [`constant`, `pcsaft`]
            - `solubility`: solubility when `constant` solubilityModel is used (mol/m^3/Pa)
    - `air` :
        - `molarMass`: (kg/mol)
    - `polymer`:
        - `heatCapacity`: (J/kg/K)
        - `polymerDensityModel`: [`constant`, `nanotools`, `pcsaft`],
        - `density`: when `constant` polymerDensityModel is used (kg/m^3)
        - `viscosityModel`: [`CastroMacosko`],
        - `maxViscosity`: maximum viscosity - for gel point detection (Pa s),
        - `molarMassNCO`: (kg/kmol)
    - `surfaceTensionModel`: [`constant` `pcsaft`],
    - `surfaceTension`: surface tension when `constant` is used (N/m)
    - `ModenaFoamViscosityModel`: use TUE model to calculate foam viscosity (t/f)
- `bubbleGrowth`:
    - `geometry`: [`3D`,`2D`]
    - `integrator`: [`dlsode`, `dlsodes`, 'cvode'] dlsodes is recommended
    - `method`: [`nonstiff`, `stiff`] stiff is recommended
    - `inertialTerm`: [t/f] use inertial term in bubble growth model
    - `solubilityCorrection`: [t/f] use solubility correction on surface tension
    - `meshCoarseningParameter`: 1.02 is recommended
    - `internalNodes`: number of nodes in shell around the bubble, 200 is recommended
    - `initialTime`: (s)
    - `finalTime`: (s)
    - `outerTimeSteps`: number of outputs
    - `maxInnerTimeSteps`: maximum number of internal time steps of integrator
    - `relativeTolerance`: 1e-8 is recommended
    - `absoluteTolerance`: 1e-8 is recommended
- `QmomKinetics`:
   - `relativeTolerance`: 1e-6 is recommended
   - `absoluteTolerance`: 1e-6 is recommended
   - `timeStep`: (s)
   - `endTime`: (s)
   - `bubbleMode`: [`mean radius`, `two nodes`]

### Preparing the test case for OpenFOAM simulation
- `0` directory: initial and boundary conditions can be defined. Further details can be found in [OpenFOAM UserGuide](http://cfd.direct/openfoam/user-guide/)
- `system` directory: includes dictioneries to control the simulation details such as time discretization, solution methods and etc. For further details, one should consult with [OpenFOAM UserGuide](http://cfd.direct/openfoam/user-guide/)
- `constant` directory: includes the mesh and and dictionaries to define the physical properties of fluid flow. There are also additional dictionaries, specifically for PU foam.

    - `kineticsProperties`: defines the details of kinetics scheme

        - `liquidMixtureDensitySurrogate`: enable the surrogate model for mixture density [on/off],
        - `blowingAgent`: name of blowing agent [n-pentane, R-11],
        - `kineticsModel`: name of kinetics model [generic, RF-1],
        - `GellingConstants`: constants for the gelling reaction

            - `A_OH`: pre-exponential factor (m^3/mol s),
            - `E_OH`: activation energy (J/mol),
            - `initCOH`: initial concentration of polyol OH groups in the mixutre (mol/m3),
            - `initCNCO`: initial concentration of isocianate NCO groups in the mixutre (mol/m3),
            - `initCW`: initial concentration of water in the mixture (mol/m3),
            - `gellingPoint`: gelling point
        - `BlowingConstants` constants for the blowing reaction

            - `A_W`: pre-exponential factor (m^3/mol s),
            - `E_W`: activation energy (J/mol)
        - `GenericConstants`: generic constants

            - `idealGasCons`: ideal gas constant (J/mol K),
            - `rhoPolymer`: density of the liquid polymer (kg/m^3),
            - `rhoBlowingAgent`: density of the blowing agent (kg/m3),
            - `molecularMassCO2`: molecular mass of carbon dioxide, (kg/kmol),
            - `molecularMassBlowingAgent`: molecular mass of blowing agent (kg/kmol),
            - `molecularMassNCO`: molecular weight of NCO (kg/kmol),
            - `molecularMassLiquidFoam`: molecular weight of liquid mixture (kg/kmol),
            - `dissolvedCO2`: weight fraction of dissolved CO2 in the mixture (kg/kg),
            - `dxdTcons`: model constant for the blowing agent (-0.01162790697 is recommended),
            - `initBlowingAgent`: initial weight fraction of blowing agent in the liquid (kg/kg),
            - `initCO2`: initial weight fraction of CO2 in the liquid (kg/kg),
            - `surfaceTension`: surface tension
       - `EnthalpyConstants`: constant for energy equation

           - `deltaOH`: reaction heat for the gelling reaction (J/mol),
           - `deltaW`: reaction heat for the blowing reaction (J/mol),
           - `PUspecificHeat`: polyurethane specific heat (J/kg K),
           - `CO2specificHeat`: CO2 specific heat (J/kg K),
           - `BGspecificHeat`: physical blowing agent in gas phase specific heat (J/kg K),
           - `BLspecificHeat`: physical blowing agent in liquid phase specific heat (J/kg K),
           - `latentHeat`: latent heat of blowing agent (J/kg)
    - `PBEProperties`: defines the details of population balance solution

        - `PBESwitch`: enbale PBE model [on/off],
        - `PBEMethod`: solution method for PBE. In the current version `QMOM` is available.
        - `nNodes`: number of quadrature approximation nodes [2-3],
        - `bubbleGrowthSurrogateSwitch`: enable bubble growth surrogate model [on/off],
        - `bubbleGrowthMode`: use of mean bubble radius or two nodes [meanRadius, twoNodes].
    - `rheologyProperties`: define the rheology model ([further details](http://onlinelibrary.wiley.com/doi/10.1002/masy.201500108/abstract))

        - `viscosityModel`: viscosity model [constant, castro-macosko, bird-carreau]
    - `simulationMode`: defines the how long the simulation should run

        - `simulationTarget`: two different modes for the simulation: mold-filling and validation. Mold-filling would disable the solution of equations when the mould is filled and validation would continue until the defined end time.

### Preparing the test case for Wall drainage simulation
- `growthRateModel`: `fromFile` uses evolution of bubble radius from Bubble growth model [`constant`,`fromFile`],
- `growthRate`: growth rate when `constant` growth rate model is used (m/s),
- `initialConditions`:
    - `centerThickness`: initial film half thickness at the center (m),
    - `domainSize`: initial size of the domain (m),
    - `filmReduction`: parameter influencing initial shape of wall and strut (0.0 = only wall, 1.0 = only strut), 1.0 is recommended,
- `physicalProperties`:
    - `viscosityModel`: `fromFile` uses evolution of viscosity from Bubble growth model [`constant`,`fromFile`],
    - `viscosity`: viscosity when `constant` viscosity model is used (Pa s),
    - `surfaceTension`: (N/m),
    - `disjoiningPressure`: parameters of disjoining pressure, consult [Schwarz and Roy (2003)][1]
        - `N`: 4.0,
        - `M`: 3.0,
        - `hst`: 1.0e-7,
        - `C`: 0.05,
        - `B1`: 0, if `B1=0` disjoining pressure is neglected
- `integration`:
    - `initialTime`: (s),
    - `timeStep`: write interval (s),
    - `outerTimeSteps`: number of outputs,
    - `method`: `stiff` is recommended [`nonstiff`, `stiff`],
    - `internalNodes`: discretization - number of finite volumes,
    - `maxInnerTimeSteps`: maximum number of internal time steps of integrator,
    - `relativeTolerance`: 1e-8 is recommended,
    - `absoluteTolerance`: 0 is recommended
- `algebraicEquationSolver`:
    - `tolerance`: 1e-4 is recommended
- `strutFilmParameter`: parameter for detection of wall/strut 1.1 is recommended

[1]: http://dx.doi.org/10.1016/S0021-9797(03)00425-9
