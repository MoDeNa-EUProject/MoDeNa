IOdictionary kineticsProperties
(
    IOobject
    (
        "kineticsProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar AOH = 
    readScalar(kineticsProperties.subDict("GellingConstants").lookup("A_OH"));
scalar EOH = 
    readScalar(kineticsProperties.subDict("GellingConstants").lookup("E_OH"));
scalar COH_0 = 
    readScalar
    (
        kineticsProperties.subDict("GellingConstants").lookup("initCOH")
    );
scalar CNCO_0 = 
    readScalar
    (
        kineticsProperties.subDict("GellingConstants").lookup("initCNCO")
    );
scalar CW_0 = 
    readScalar
    (
        kineticsProperties.subDict("GellingConstants").lookup("initCW")
    );
scalar XOH_Gel = 
    readScalar
    (
        kineticsProperties.subDict("GellingConstants").lookup("gellingPoint")
    );

if (EOH < 0.0 || COH_0 < 0.0 || CNCO_0 < 0.0 || CW_0 < 0.0 || XOH_Gel < 0.0 )
{
    FatalErrorIn(args.executable())
    << "Invalid value for gelling constants!"
    "\nGelling constats cannot be negative."
    << exit(FatalError);
}

word blowingAgent(kineticsProperties.lookup("blowingAgent"));

if
(
    !(blowingAgent == "n-pentane" 
    || blowingAgent == "R-11" 
    || blowingAgent == "no")
)
{
    FatalErrorIn(args.executable())
    << blowingAgent << " Is an invalid type of physical blowing agent."
    << "\nValid types are: "
    << "\nn-pentane,"
    << "\nR-11"
    << "\nno"
    << exit(FatalError);
}

scalar AW = 
    readScalar(kineticsProperties.subDict("BlowingConstants").lookup("A_W"));
scalar EW = 
    readScalar(kineticsProperties.subDict("BlowingConstants").lookup("E_W"));

if (EW < 0.0)
{
    FatalErrorIn(args.executable())
    << "Invalid value for blowing constants!"
    "\nActivation energy cannot be negative."
    << exit(FatalError);
}

scalar RR = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("idealGasCons")
    );
scalar rhoPoly = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("rhoPolymer")
    );
scalar rhoBL = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("rhoBlowingAgent")
    );
scalar M_CO2 = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("molecularMassCO2")
    );
scalar M_B = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("molecularMassBlowingAgent")
    );
scalar M_NCO = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("molecularMassNCO")
    );
scalar CO2_D = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("dissolvedCO2")
    );
scalar dxdT = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("dxdTcons")
    );
scalar L0 = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("initBlowingAgent")
    );
scalar CO2_0 = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("initCO2")
    );
scalar M_liq = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("molecularMassLiquidFoam")
    );
scalar surfaceTension = 
    readScalar
    (
        kineticsProperties.subDict("GenericConstants").lookup("surfaceTension")
    );

if
(
    RR < 0.0 || rhoPoly < 0.0 || rhoBL < 0.0 
    || M_B < 0.0 || M_CO2 < 0.0 || M_NCO < 0.0 
    || L0 < 0.0 || M_liq < 0.0 || surfaceTension < 0.0
)
{
    FatalErrorIn(args.executable())
    << "\nInvalid value for generic constants!"
    << exit(FatalError);
}

// Enthalpy
scalar DH_OH = 
    readScalar
    (
        kineticsProperties.subDict("EnthalpyConstants").lookup("deltaOH")
    );
scalar DH_W = 
    readScalar
    (
        kineticsProperties.subDict("EnthalpyConstants").lookup("deltaW")
    );
scalar C_Poly = 
    readScalar
    ( 
        kineticsProperties.subDict("EnthalpyConstants").lookup("PUspecificHeat")
    );
scalar C_CO2 = 
    readScalar
    (
        kineticsProperties.subDict("EnthalpyConstants").lookup("CO2specificHeat")
    );
scalar C_BG = 
    readScalar
    (
        kineticsProperties.subDict("EnthalpyConstants").lookup("BGspecificHeat")
    );
scalar C_BL = 
    readScalar
    (
        kineticsProperties.subDict("EnthalpyConstants").lookup("BLspecificHeat")
    );
scalar latenth = 
    readScalar
    ( 
        kineticsProperties.subDict("EnthalpyConstants").lookup("latentHeat")
    );
// switch for the liquid mixture density
    Switch liquidMixtureDensitySurrogate
    (
        kineticsProperties.lookup("liquidMixtureDensitySurrogate")
    );
Info<< "Reading field XW \n" << endl;
volScalarField XW
(
    IOobject
    (
        "XW",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);

Info<< "Creating field rho_foam*alpha2\n" << endl;
volScalarField rhoFoamAlpha2 
(
    "rhoFoamAlpha2",
    rho_foam*alpha2
);

Info<< "Creating surface field of rho_foam*alpha2*U\n" << endl;
surfaceScalarField frhoFoamAlpha2
(
    fvc::interpolate(rhoFoamAlpha2)*phi
);

volScalarField::DimensionedInternalField XWSource
(
    IOobject
    (
        "XWSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("XWSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field XOH \n" << endl;
volScalarField XOH
(
    IOobject
    (
        "XOH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);

Info<< "\nCreating source fields ScNonXOH for XOH: " << endl;
volScalarField::DimensionedInternalField ScNonXOH
(
    IOobject
    (
        "ScNonXOH",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("ScNonXOH", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "\nCreating first order coefficient for XOH source: " << endl;
volScalarField::DimensionedInternalField Coeff_1stOrd
(
    IOobject
    (
        "Coeff_1stOrd",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("Coeff_1stOrd", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Creating SuXOH and SpXOH for the second order term of in XOH source" << endl;
volScalarField::DimensionedInternalField SuXOH
(
    IOobject
    (
        "SuXOH",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("SuXOH", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);
volScalarField::DimensionedInternalField SpXOH
(
    IOobject
    (
        "SpXOH",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("SpXOH", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

volScalarField::DimensionedInternalField TSource 
(
    IOobject
    (
        "TSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("TSource", dimensionSet(1,-3,-1,1,0,0,0), 0.0)
);

Info<< "Reading field wBA_g, blowing agent in gas \n" << endl;
volScalarField wBA_g
(
    IOobject
    (
       "wBA_g",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);
volScalarField::DimensionedInternalField wBA_gSource
(
    IOobject
    (
        "wBA_gSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("wBA_gSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field wBA_l, blowing agent in liquid \n" << endl;
volScalarField wBA_l
(
    IOobject
    (
        "wBA_l",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);

volScalarField::DimensionedInternalField wBA_lSource
(
    IOobject
    (
        "wBA_lSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("wBA_lSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field wCO2_l, CO2 in liquid \n" << endl;
volScalarField wCO2_l
(
    IOobject
    (
        "wCO2_l",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);
volScalarField::DimensionedInternalField wCO2_lSource
(
    IOobject
    (
        "wCO2_lSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("wCO2_lSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field wCO2_g, CO2 in gas \n" << endl;
volScalarField wCO2_g
(
    IOobject
    (
        "wCO2_g",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);
volScalarField::DimensionedInternalField wCO2_gSource
(
    IOobject
    (
        "wCO2_gSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("wCO2_gSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field g1_CO2, source of mOne due to the generation of CO2\n" << endl;
volScalarField g1_CO2
(
    IOobject
    (
        "g1_CO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);
Info<< "Reading field g1_BA, source of mOne due to the evaporation of blowing agent\n" << endl;
volScalarField g1_BA
(
    IOobject
    (
        "g1_BA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh    
);
Info<< "Reading field cc1, a constant based on the growth rate due to blowing agent\n" << endl;
volScalarField cc1
(
    IOobject
    (
        "cc1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh     
);

volScalarField CO2
(
    IOobject
    (
        "CO2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CO2", dimless, 0.0)
);
volScalarField creamT
(
    IOobject
    (
        "creamT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("creamT", dimless, 0.0)
);
volScalarField TS
(
    IOobject
    (
        "TS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TS", dimensionSet(0,0,0,1,0,0,0), 298.0)
);
volScalarField::DimensionedInternalField TSSource 
(
    IOobject
    (
        "TSSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("TSSource", dimensionSet(0,0,-1,1,0,0,0), 0.0)
);