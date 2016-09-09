IOdictionary PBEProperties
(
    IOobject
    (
        "PBEProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch PBESwitch
(
    PBEProperties.lookup("PBESwitch")
);
Switch bubbleGrowthSurrogateSwitch
    (
        PBEProperties.lookup("bubbleGrowthSurrogateSwitch")
    );
word PBEMethod(PBEProperties.lookup("PBEMethod"));
if( !(PBEMethod == "QMOM") )
{
    FatalErrorIn(args.executable())
    << "Invalid PBE solution method " << PBEMethod
    << exit(FatalError);
}

word bubbleGrowthMode(PBEProperties.lookup("bubbleGrowthMode"));
if
(
    !( bubbleGrowthMode == "twoNodes"
    || bubbleGrowthMode == "meanRadius"
    )
)
{
    FatalErrorIn(args.executable())
    << "Invalid bubble growth mode " << bubbleGrowthMode
    << "\nValid modes are: "
    << "\ntwoNodes,"
    << "\nmeanRadius"
    << exit(FatalError);
}

scalar nNodes = readScalar(PBEProperties.lookup("nNodes"));
if
(
    !(nNodes == 2
    || nNodes == 3)
)
{
    FatalErrorIn(args.executable())
    << "Invalid number of quadrature nodes " << nNodes
    << exit(FatalError);
}

volScalarField cellvolume
(
    IOobject
    (
        "cellvolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimVolume, 0.0)
);

volVectorField cellcentroid
(
    IOobject
    (
        "cellcentroid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimLength, vector::zero)
);

Info<< "Reading field mZero \n" << endl;
volScalarField mZero
(
    IOobject
    (
        "mZero",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ,
         IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField::Internal mZeroSource
(
    IOobject
    (
        "mZeroSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mZeroSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

surfaceScalarField phic(mag(phi/mesh.magSf()));
phic = min(interface.cAlpha()*phic, max(phic));
surfaceScalarField phir(phic*interface.nHatf());

surfaceScalarField alpha1phir
(
    fvc::interpolate(alpha1)*phir
);

Info<< "Reading field mOne \n" << endl;
volScalarField mOne
(
    IOobject
    (
        "mOne",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField::Internal mOneSource
(
    IOobject
    (
        "mOneSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mOneSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field mTwo \n" << endl;
volScalarField mTwo
(
    IOobject
    (
        "mTwo",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField::Internal mTwoSource
(
    IOobject
    (
        "mTwoSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mTwoSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field mThree \n" << endl;
volScalarField mThree
(
    IOobject
    (
        "mThree",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField::Internal mThreeSource
(
    IOobject
    (
        "mThreeSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mThreeSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field mFour \n" << endl;
volScalarField mFour
(
    IOobject
    (
        "mFour",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField::Internal mFourSource
(
    IOobject
    (
        "mFourSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mFourSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field mFive \n" << endl;
volScalarField mFive
(
   IOobject
   (
        "mFive",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
   ),
   mesh
);

volScalarField::Internal mFiveSource
(
    IOobject
    (
        "mFiveSource",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("mFiveSource", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field M0 \n" << endl;
volScalarField M0
(
    IOobject
    (
        "M0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field M1 \n" << endl;
volScalarField M1
(
    IOobject
    (
        "M1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field M2 \n" << endl;
volScalarField M2
(
    IOobject
    (
        "M2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field M3 \n" << endl;
volScalarField M3
(
    IOobject
    (
        "M3",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field M4 \n" << endl;
volScalarField M4
(
    IOobject
    (
        "M4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field M5 \n" << endl;
volScalarField M5
(
    IOobject
    (
        "M5",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fields weights and nodes \n" << endl;
volScalarField weight0
(
    IOobject
    (
        "weight0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh

);

volScalarField weight1
(
    IOobject
    (
        "weight1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField node0
(
    IOobject
    (
        "node0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField node1
(
    IOobject
    (
        "node1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField weight2
(
    IOobject
    (
        "weight2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField node2
(
    IOobject
    (
        "node2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "End of \"MomFields\"\n" << endl;
// End of MomFields