{
    dimensionedScalar continuityError;
    continuityError =
        mag(initialFoamMass - fvc::domainIntegrate(rho_foam*alpha2))
        /initialFoamMass*100.0;
    cumulativeContinuityError = cumulativeContinuityError + continuityError;

    Info<< "Continuity error: "
        << continuityError.value() << " %" << endl;
}