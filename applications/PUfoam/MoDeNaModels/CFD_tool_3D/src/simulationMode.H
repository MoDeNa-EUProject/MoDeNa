IOdictionary simulationMode
(
    IOobject
    (
        "simulationMode",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word simulationTarget(simulationMode.lookup("simulationTarget"));

if
(
    !(simulationTarget == "mold-filling"
   || simulationTarget == "validation")
)
{
    FatalErrorIn(args.executable())
    << "Invalid simulation target " << simulationTarget
    << "\nValid simulation targets are: "
    << "\nmold-filling,"
    << "\nvalidation,"
    << exit(FatalError);
}