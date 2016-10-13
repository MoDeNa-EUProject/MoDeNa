## postProcessBSD
This is a post-processing application to analyze the results of population balance equation. It calcualtes the mean and variance of the bubble size distribution at each time step and write the results as field variables in each time directory. They can be further post-processed using `postProcess -func probe` command to probe the values on the target locations of the domain.

### Dependencies and Compilation
In order to run this application OpenFOAM 4.0 should be installed. The installation of OpenFOAM 4.0 has been described in `foamExpansion` directory. The compilation process for this application includes the following:

```
cd $HOME/MoDeNa/applications/PUfoam/MoDeNaModels/CFD_tool_3D/applications/ postProcessBSD
wclean
wmake
```

### How to run?
After running the `workflow_3D` in `foamExpansion` directory:

```
cd #latest launcher directory where the CFD results are stored
cp -r $FOAM_ETC/caseDicts/postProcessing/probes/probes system
cd system
```
Now open the `probes` dictionary with an editor to specify the required fields (e.g. BSDMean, BSDVariance) at target points. Here is an example of `probes` dictionary for probing `BSDMean` and `BSDVariance`:
```
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Web:      www.OpenFOAM.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Writes out values of fields from cells nearest to specified locations.

\*---------------------------------------------------------------------------*/

#includeEtc "caseDicts/postProcessing/probes/probes.cfg"

fields (
        BSDMean
        BSDVariance
       );
probeLocations
(
    (0.05 0.01 0.005)
);

// ************************************************************************* //
```
This process will be complete by adding the `#includeFunc` directive to `functions` in the `controlDict` file in `system` directory:

```
functions
{
    #includeFunc  probes

}
```

Then you can run `postProcessBSD` application and probe the results by typing:
```
postProcessBSD
postProcess -func probes
```
