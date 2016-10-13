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
   || KineticsModel == "RF-1" )
)
    {
        FatalErrorIn(args.executable())
            << "Invalid kinetics model " << KineticsModel
            << "\nValid kinetics models are: "
            << "\ngeneric,"
            << "\nRF-1,"
            << exit(FatalError);
    }
scalar AOH, EOH, COH_0, CNCO_0, CW_0, AW, EW;
AOH=EOH=COH_0=CNCO_0=CW_0=AW=EW=0.0;
if (KineticsModel == "generic")
{
    AOH =
        readScalar(kineticsProperties.subDict("GellingConstants").lookup("A_OH"));
    EOH =
        readScalar(kineticsProperties.subDict("GellingConstants").lookup("E_OH"));

    COH_0 =
        readScalar
        (
            kineticsProperties.subDict("GellingConstants").lookup("initCOH")
        );
    CNCO_0 =
        readScalar
        (
            kineticsProperties.subDict("GellingConstants").lookup("initCNCO")
        );
    CW_0 =
        readScalar
        (
            kineticsProperties.subDict("GellingConstants").lookup("initCW")
        );
    AW =
        readScalar(kineticsProperties.subDict("BlowingConstants").lookup("A_W"));
    EW =
        readScalar(kineticsProperties.subDict("BlowingConstants").lookup("E_W"));
    if (EOH < 0.0 || COH_0 < 0.0 || CNCO_0 < 0.0 || CW_0 < 0.0 )
    {
        FatalErrorIn(args.executable())
        << "Invalid value for gelling constants!"
        "\nGelling constats cannot be negative."
        << exit(FatalError);
    }
    if (EW < 0.0)
    {
        FatalErrorIn(args.executable())
        << "Invalid value for blowing constants!"
        "\nActivation energy cannot be negative."
        << exit(FatalError);
    }
}

scalar XOH_Gel =
    readScalar
    (
        kineticsProperties.subDict("GellingConstants").lookup("gellingPoint")
    );
if (XOH_Gel < 0.0)
{
    FatalErrorIn(args.executable())
    << "Invalid value for gelling point!"
    "\nGelling point cannot be negative."
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

volScalarField::Internal XWSource
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

volScalarField::Internal XOHSource
(
    IOobject
    (
        "XOHSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("XOHSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
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
volScalarField::Internal ScNonXOH
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
volScalarField::Internal Coeff_1stOrd
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
volScalarField::Internal SuXOH
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
volScalarField::Internal SpXOH
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

volScalarField::Internal TSource
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
volScalarField::Internal wBA_gSource
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

volScalarField::Internal wBA_lSource
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
volScalarField::Internal wCO2_lSource
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
volScalarField::Internal wCO2_gSource
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
volScalarField::Internal TSSource
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
volScalarField source_Catalyst_1
    (
         IOobject
        (
         "source_Catalyst_1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_Catalyst_1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_A0" << endl;
volScalarField source_CE_A0
    (
         IOobject
        (
         "source_CE_A0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_A0", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_A1" << endl;
volScalarField source_CE_A1
    (
         IOobject
        (
         "source_CE_A1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_A1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_B" << endl;
volScalarField source_CE_B
    (
         IOobject
        (
         "source_CE_B",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_B", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_B2" << endl;
volScalarField source_CE_B2
    (
         IOobject
        (
         "source_CE_B2",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_B2", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_I0" << endl;
volScalarField source_CE_I0
    (
         IOobject
        (
         "source_CE_I0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_I0", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_I1" << endl;
volScalarField source_CE_I1
    (
         IOobject
        (
         "source_CE_I1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_I1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_I2" << endl;
volScalarField source_CE_I2
    (
         IOobject
        (
         "source_CE_I2",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_I2", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_PBA" << endl;
volScalarField source_CE_PBA
    (
         IOobject
        (
         "source_CE_PBA",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_PBA", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Breac" << endl;
volScalarField source_CE_Breac
    (
         IOobject
        (
         "source_CE_Breac",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Breac", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Areac0" << endl;
volScalarField source_CE_Areac0
    (
         IOobject
        (
         "source_CE_Areac0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Areac0", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Areac1" << endl;
volScalarField source_CE_Areac1
    (
         IOobject
        (
         "source_CE_Areac1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Areac1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Ireac0" << endl;
volScalarField source_CE_Ireac0
    (
         IOobject
        (
         "source_CE_Ireac0",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Ireac0", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Ireac1" << endl;
volScalarField source_CE_Ireac1
    (
         IOobject
        (
         "source_CE_Ireac1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Ireac1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_CE_Ireac2" << endl;
volScalarField source_CE_Ireac2
    (
         IOobject
        (
         "source_CE_Ireac2",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_CE_Ireac2", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_Bulk" << endl;
volScalarField source_Bulk
    (
         IOobject
        (
         "source_Bulk",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_Bulk", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1" << endl;
volScalarField source_R_1
    (
         IOobject
        (
         "source_R_1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_R_1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_mass" << endl;
volScalarField source_R_1_mass
    (
         IOobject
        (
         "source_R_1_mass",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("source_R_1_mass", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_temp_RF1" << endl;
volScalarField source_R_1_temp_RF1
    (
         IOobject
        (
         "source_R_1_temp_RF1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("sourceR_1_temp_RF1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
Info<< "Creating source term source_R_1_vol_RF1" << endl;
volScalarField source_R_1_vol_RF1
    (
         IOobject
        (
         "source_R_1_vol_RF1",
         runTime.timeName(),
         mesh
        ),
        mesh,
        dimensionedScalar("sourceR_1_vol_RF1", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

if (KineticsModel == "RF-1")
{
    forAll(mesh.C(), celli)
    {
        if (alpha1[celli] == 0.0)
        {
            scalar foam_index;
            foam_index = celli;
            CW_0 = scalar(1000.0)*CE_B2[foam_index];
            COH_0 = scalar(1000.0)*(CE_A0[foam_index] + CE_A1[foam_index]);
            break;
        }
    }
    reduce(CW_0,maxOp<scalar>());
    reduce(COH_0,maxOp<scalar>());
}
