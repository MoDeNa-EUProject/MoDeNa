/**
* @ingroup mod_3Dcfd
* @file rheologyFields.H
* @brief creates the field variables for the rheology
*
* @detail
* This file creates the required variables for the rheology model. The
* `rheologyProperties` dictionary is where the user can define the
* `viscosityModel`. Four choices has been identified. The first one uses
* a constant value and treats the PU foam similar to water. The
* `castro-macosko` model applies a Newtonian approach relating the foam
* viscosity to the temperature and conversion. The thirs choice `bird-carreau`
* ,however , computes the foam apparent viscocity based on a non-Newtonian
* method. The last option, `MoDeNaRheology`, applies a surrogate model for the
* foam apparent viscosity.
*/

///@cond
IOdictionary rheologyProperties
(
    IOobject
    (
        "rheologyProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word viscosityModel(rheologyProperties.lookup("viscosityModel"));

if
(
    !(viscosityModel == "constant"
   || viscosityModel == "castro-macosko"
   || viscosityModel == "bird-carreau"
   || viscosityModel == "MoDeNaRheology")
)
{
    FatalErrorIn(args.executable())
    << "Invalid viscosity model " << viscosityModel
    << "\nValid viscosity models are: "
    << "\nconstant,"
    << "\ncastro-macosko,"
    << "\nbird-carreau,"
    << "\nMoDeNaRheology"
    << exit(FatalError);
}

volScalarField muFoam
(
    IOobject
    (
        "muFoam",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muFoam", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
   );
volScalarField muAir
(
    IOobject
    (
        "muAir",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muAir", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
   );
volScalarField muMixture
(
    IOobject
    (
        "muMixture",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muMixture", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
   );

volScalarField XNCO
(
    IOobject
    (
        "XNCO",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volScalarField mu0
(
    IOobject
    (
        "mu0",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mu0", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
);
volScalarField muinf
(
    IOobject
    (
        "muinf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muinf", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
);
volScalarField muFoamCorr
(
    IOobject
    (
        "muFoamCorr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("muFoamCorr", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
);
///@endcond