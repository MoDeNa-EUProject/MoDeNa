    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
    twoPhaseMixtureThermo mixture(mesh);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    Info<< "Reading thermophysical properties\n" << endl;

    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();

    Info<< "Reading field rho_foam \n" << endl;
    volScalarField rho_foam
    (
        IOobject
        (
	    "rho_foam",
	     runTime.timeName(),
	     mesh,
	     IOobject::MUST_READ,
	     IOobject::AUTO_WRITE
        ),
         mesh
    );

    Info<< "Reading field rho_gas \n" << endl;
    volScalarField rho_gas
    (
        IOobject
        (
	    "rho_gas",
	     runTime.timeName(),
	     mesh,
	     IOobject::MUST_READ,
	     IOobject::AUTO_WRITE
        ),
         mesh
    );
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho_gas + (1-alpha1)*rho_foam
    );
    volScalarField Psi1
    (
        IOobject
        (
	    "Psi1",
	     runTime.timeName(),
	     mesh,
	     IOobject::READ_IF_PRESENT,
	     IOobject::AUTO_WRITE
        ),
         mesh,
	    dimensionedScalar("Psi1", dimensionSet(0,-2,2,0,0,0,0), 1e-5)
    );

    volScalarField Psi2
    (
        IOobject
        (
	    "Psi2",
	     runTime.timeName(),
	     mesh,
	     IOobject::READ_IF_PRESENT,
	     IOobject::AUTO_WRITE
        ),
         mesh,
	    dimensionedScalar("Psi2", dimensionSet(0,-2,2,0,0,0,0), 1e-2)
    );

    volScalarField Psi_bubble
    (
        IOobject
        (
        "Psi_bubble",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
        ),
         mesh,
        dimensionedScalar("Psi_bubble", dimensionSet(0,-2,2,0,0,0,0), 1e-2)
    );

    volScalarField Psi_foam
    (
        IOobject
        (
        "Psi_foam",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
        ),
         mesh,
        dimensionedScalar("Psi_foam", dimensionSet(0,-2,2,0,0,0,0), 1e-2)
    );

    dimensionedScalar pMin(mixture.lookup("pMin"));

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );

    volScalarField dgdt
    (
        pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
    );


    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, mixture);
    #include "MomFields.H"
    #include "KineticsFields.H"
    #include "rheologyFields.H"
    #include "thermalFields.H"
    #include "simulationMode.H"
    // Construct compressible turbulence model
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New(rho, U, rhoPhi, mixture)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));

    dimensionedScalar densityDimension ("densityDimension", dimensionSet(1,-3,0,0,0,0,0), 1.0);
    dimensionedScalar initialFoamMass =
        (
            fvc::domainIntegrate(rho_foam*alpha2)
            + densityDimension*L0*rhoPoly*fvc::domainIntegrate(alpha2)
        );
    dimensionedScalar cumulativeContinuityError = 0.0;