/**
* @ingroup mod_3Dcfd
* @file continuityErrors.H
* @brief calculates and reports the continuity error
*
* @detail
* This file calculates the continuity error based on the initial mass of foam
* `initialFoamMass`  and the amount of foam at each time  step`fvc::domainIntegrate(rho_foam*alpha2)`.
*/
///@cond
{
    dimensionedScalar continuityError;
    continuityError =
        mag(initialFoamMass - fvc::domainIntegrate(rho_foam*alpha2))
        /initialFoamMass*100.0;
    cumulativeContinuityError = cumulativeContinuityError + continuityError;

    Info<< "Continuity error: "
        << continuityError.value() << " %" << endl;
}
///@endcond