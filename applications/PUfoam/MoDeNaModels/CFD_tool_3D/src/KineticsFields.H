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
word KineticsModel(kineticsProperties.lookup("KineticsModel"));

if 
( 
    !(KineticsModel == "generic" 
   || KineticsModel == "simple" 
   || KineticsModel == "RF-1" ) 
)
    {
        FatalErrorIn(args.executable())
            << "Invalid kinetics model " << KineticsModel
            << "\nValid kinetics models are: "
            << "\ngeneric,"
            << "\nsimple,"
            << "\nRF-1,"
            << exit(FatalError);
    }

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
// Required fields for generic kinetics model
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
// Required fields for simple kinetics surrogate models
Info<< "Reading field EG_NCO \n" << endl;
    volScalarField EG_NCO
    (
        IOobject
        (
            "EG_NCO",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field EG_OH \n" << endl;
    volScalarField EG_OH
    (
        IOobject
        (
            "EG_OH",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field H2O \n" << endl;
    volScalarField H2O
    (
        IOobject
        (
            "H2O",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field CO2_Sur \n" << endl;
    volScalarField CO2_Sur
    (
        IOobject
        (
            "CO2_Sur",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field PENTANE \n" << endl;
    volScalarField PENTANE
    (
        IOobject
        (
            "PENTANE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field POLYMER \n" << endl;
    volScalarField POLYMER
    (
        IOobject
        (
            "POLYMER",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field POLMERBLOW \n" << endl;
    volScalarField POLMERBLOW
    (
        IOobject
        (
            "POLMERBLOW",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field UREA \n" << endl;
    volScalarField UREA
    (
        IOobject
        (
            "UREA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );

Info<< "Reading field R_1_temp \n" << endl;
    volScalarField R_1_temp
    (
        IOobject
        (
            "R_1_temp",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
         mesh     
    );

// Source terms for the simple kinetics surrogate model
Info<< "Creating source term source_EG_NCO" << endl;
volScalarField::DimensionedInternalField source_EG_NCO
    (
         IOobject
        (
            "source_EG_NCO",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_EG_NCO", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_EG_OH" << endl;
volScalarField::DimensionedInternalField source_EG_OH
    (
         IOobject
        (
            "source_EG_OH",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_EG_OH", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_H2O" << endl;
volScalarField::DimensionedInternalField source_H2O
    (
         IOobject
        (
            "source_H2O",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_H2O", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CO2" << endl;
volScalarField::DimensionedInternalField source_CO2
    (
         IOobject
        (
            "source_CO2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_CO2", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_PENTANE" << endl;
volScalarField::DimensionedInternalField source_PENTANE
    (
         IOobject
        (
            "source_PENTANE",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_PENTANE", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_POLYMER" << endl;
volScalarField::DimensionedInternalField source_POLYMER
    (
         IOobject
        (
            "source_POLYMER",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_POLYMER", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_POLMERBLOW" << endl;
volScalarField::DimensionedInternalField source_POLMERBLOW
    (
         IOobject
        (
            "source_POLMERBLOW",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_POLMERBLOW", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_UREA" << endl;
volScalarField::DimensionedInternalField source_UREA
    (
         IOobject
        (
            "source_UREA",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_UREA", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_temp" << endl;
volScalarField::DimensionedInternalField source_R_1_temp
    (
         IOobject
        (
            "source_R_1_temp",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("source_R_1_temp", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
// Required fields for RF-1 kinetics surrogate model
Info<< "Reading field Catalyst_1 \n" << endl;
    volScalarField Catalyst_1
    (
        IOobject
        (
            "Catalyst_1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_A0 \n" << endl;
    volScalarField CE_A0
    (
        IOobject
        (
            "CE_A0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_A1 \n" << endl;
    volScalarField CE_A1
    (
        IOobject
        (
            "CE_A1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_B \n" << endl;
    volScalarField CE_B
    (
        IOobject
        (
            "CE_B",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_B2 \n" << endl;
    volScalarField CE_B2
    (
        IOobject
        (
            "CE_B2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_I0 \n" << endl;
    volScalarField CE_I0
    (
        IOobject
        (
            "CE_I0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_I1 \n" << endl;
    volScalarField CE_I1
    (
        IOobject
        (
            "CE_I1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_I2 \n" << endl;
    volScalarField CE_I2
    (
        IOobject
        (
            "CE_I2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_PBA \n" << endl;
    volScalarField CE_PBA
    (
        IOobject
        (
            "CE_PBA",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_Breac \n" << endl;
    volScalarField CE_Breac
    (
        IOobject
        (
            "CE_Breac",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_Areac0 \n" << endl;
    volScalarField CE_Areac0
    (
        IOobject
        (
            "CE_Areac0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_Areac1 \n" << endl;
    volScalarField CE_Areac1
    (
        IOobject
        (
            "CE_Areac1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_Ireac0 \n" << endl;
    volScalarField CE_Ireac0
    (
        IOobject
        (
            "CE_Ireac0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_Ireac1 \n" << endl;
    volScalarField CE_Ireac1
    (
        IOobject
        (
            "CE_Ireac1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field CE_Ireac2 \n" << endl;
    volScalarField CE_Ireac2
    (
        IOobject
        (
            "CE_Ireac2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field Bulk \n" << endl;
    volScalarField Bulk
    (
        IOobject
        (
            "Bulk",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field R_1 \n" << endl;
    volScalarField R_1
    (
        IOobject
        (
            "R_1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field R_1_mass \n" << endl;
    volScalarField R_1_mass
    (
        IOobject
        (
            "R_1_mass",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field R_1_temp_RF1 \n" << endl;
    volScalarField R_1_temp_RF1
    (
        IOobject
        (
            "R_1_temp_RF1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
Info<< "Reading field R_1_vol_RF1 \n" << endl;
    volScalarField R_1_vol_RF1
    (
        IOobject
        (
            "R_1_vol_RF1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh     
    );
// source terms for the RF-1 kinetics surrogate models
Info<< "Creating source term source_Catalyst_1" << endl;
volScalarField::DimensionedInternalField source_Catalyst_1
    (
         IOobject
        (
         "source_Catalyst_1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_Catalyst_1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_A0" << endl;
volScalarField::DimensionedInternalField source_CE_A0
    (
         IOobject
        (
         "source_CE_A0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_A0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_A1" << endl;
volScalarField::DimensionedInternalField source_CE_A1
    (
         IOobject
        (
         "source_CE_A1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_A1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_B" << endl;
volScalarField::DimensionedInternalField source_CE_B
    (
         IOobject
        (
         "source_CE_B",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_B", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_B2" << endl;
volScalarField::DimensionedInternalField source_CE_B2
    (
         IOobject
        (
         "source_CE_B2",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_B2", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_I0" << endl;
volScalarField::DimensionedInternalField source_CE_I0
    (
         IOobject
        (
         "source_CE_I0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_I0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_I1" << endl;
volScalarField::DimensionedInternalField source_CE_I1
    (
         IOobject
        (
         "source_CE_I1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_I1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_I2" << endl;
volScalarField::DimensionedInternalField source_CE_I2
    (
         IOobject
        (
         "source_CE_I2",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_I2", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_PBA" << endl;
volScalarField::DimensionedInternalField source_CE_PBA
    (
         IOobject
        (
         "source_CE_PBA",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_PBA", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Breac" << endl;
volScalarField::DimensionedInternalField source_CE_Breac
    (
         IOobject
        (
         "source_CE_Breac",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Breac", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Areac0" << endl;
volScalarField::DimensionedInternalField source_CE_Areac0
    (
         IOobject
        (
         "source_CE_Areac0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Areac0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Areac1" << endl;
volScalarField::DimensionedInternalField source_CE_Areac1
    (
         IOobject
        (
         "source_CE_Areac1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Areac1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Ireac0" << endl;
volScalarField::DimensionedInternalField source_CE_Ireac0
    (
         IOobject
        (
         "source_CE_Ireac0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Ireac0", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Ireac1" << endl;
volScalarField::DimensionedInternalField source_CE_Ireac1
    (
         IOobject
        (
         "source_CE_Ireac1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Ireac1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Ireac2" << endl;
volScalarField::DimensionedInternalField source_CE_Ireac2
    (
         IOobject
        (
         "source_CE_Ireac2",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Ireac2", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_Bulk" << endl;
volScalarField::DimensionedInternalField source_Bulk
    (
         IOobject
        (
         "source_Bulk",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_Bulk", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1" << endl;
volScalarField::DimensionedInternalField source_R_1
    (
         IOobject
        (
         "source_R_1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_R_1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_mass" << endl;
volScalarField::DimensionedInternalField source_R_1_mass
    (
         IOobject
        (
         "source_R_1_mass",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_R_1_mass", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_temp_RF1" << endl;
volScalarField::DimensionedInternalField source_R_1_temp_RF1
    (
         IOobject
        (
         "source_R_1_temp_RF1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("sourceR_1_temp_RF1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_vol_RF1" << endl;
volScalarField::DimensionedInternalField source_R_1_vol_RF1
    (
         IOobject
        (
         "source_R_1_vol_RF1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("sourceR_1_vol_RF1", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    );