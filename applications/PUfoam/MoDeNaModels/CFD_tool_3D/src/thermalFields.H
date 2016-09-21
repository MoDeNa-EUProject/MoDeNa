volScalarField thermalConductivity
(
    IOobject
    (
        "thermalConductivity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalConductivity", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.024)
   );

Switch thermalConductivitySurrogateSwitch
    (
        kineticsProperties.lookup("thermalConductivitySurrogateSwitch")
    );

volScalarField::Internal XNCOSource
(
    IOobject
    (
        "XNCOSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("XNCOSource", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
);

volScalarField thermalDiff_foam
(
    IOobject
    (
        "thermalDiff_foam",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDiff_foam", dimensionSet(0, 2, -1, 0, 0, 0, 0), 4e-5)
);
volScalarField thermalDiffusivity
(
    IOobject
    (
        "thermalDiffusivity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDiffusivity", dimensionSet(0, 2, -1, 0, 0, 0, 0), 4e-5)
);
dimensionedScalar kappadim
(
    "kappadim",
    dimensionSet (1, 1, -3, -1, 0, 0, 0),
    1
);
dimensionedScalar cpAir
(
    "cpAir",
    dimensionSet (0, 2, -2, -1, 0, 0, 0),
    4186
);

volScalarField thermalDiff_gas
(
    IOobject
    (
        "thermalDiff_gas",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDiff_gas", dimensionSet(0, 2, -1, 0, 0, 0, 0), 1e-1)
);