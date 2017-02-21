/**
* @ingroup mod_3Dcfd
* @file simulationMode.H
* @brief defines mold-filling and validation as two simulation modes
*
* @detail
* It has been foreseen that the simulation of polyurethane foam can be run
* in two different modes including the __mold-filling__ and __validation__.
* The first mode assumes the solution of equations stops when the
* computational domain is filled, while the second option allows user to run
* the simulation for a targeted time.
* This can be defined with a dictionary named `simulationMode`.
*
*/

///@cond
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
///@endcond