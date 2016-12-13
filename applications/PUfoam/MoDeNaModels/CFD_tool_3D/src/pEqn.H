{
    Psi1 = rho_gas/p_rgh;
    Psi2 = rho_foam/p_rgh;

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
            //fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    if (pimple.transonic())
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(Psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(Psi2)*phi);

        p_rghEqnComp1 =
            fvc::ddt(rho_gas) + fvc::div(phi, rho_gas)
          - fvc::Sp(fvc::div(phi), rho_gas)
          + correction
            (
                Psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp1.ref().faceFluxCorrectionPtr());
        p_rghEqnComp1.ref().relax();

        p_rghEqnComp2 =
            fvc::ddt(rho_foam) + fvc::div(phi, rho_foam)
          - fvc::Sp(fvc::div(phi), rho_foam)
          + correction
            (
                Psi2*fvm::ddt(p_rgh)
              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp2.ref().faceFluxCorrectionPtr());
        p_rghEqnComp2.ref().relax();
    }

    else
    {
        p_rghEqnComp1 =
        (
            fvc::ddt(alpha1, rho_gas)
          + fvc::div(fvc::interpolate(alpha1)*phi, rho_gas)
          - fvc::Sp(fvc::ddt(alpha1)
          + fvc::div(fvc::interpolate(alpha1)*phi), rho_gas)
        )/rho_gas
      + (alpha1*Psi1/rho_gas)*correction(fvm::ddt(p_rgh));

        p_rghEqnComp2 =
        (
            fvc::ddt(alpha2, rho_foam)
          + fvc::div(fvc::interpolate(alpha2)*phi, rho_foam)
          - fvc::Sp(fvc::ddt(alpha2)
          + fvc::div(fvc::interpolate(alpha2)*phi), rho_foam)
        )/rho_foam
      + (alpha2*Psi2/rho_foam)*correction(fvm::ddt(p_rgh));
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );
        solve
        (
            (
                p_rghEqnComp1() + p_rghEqnComp2()
            )
          + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );
        if (pimple.finalNonOrthogonalIter())
        {
            p = max(p_rgh + (alpha1*rho_gas + (1.0 - alpha1)*rho_foam)*gh, pMin);
            p_rgh = p - (alpha1*rho_gas + (1.0 - alpha1)*rho_foam)*gh;
            dimensionedScalar dumdgdt ("dumdgdt",dimensionSet(1,-3,0,0,0,0,0), 1.0);
            dgdt =
            dumdgdt*
            (
                pos((1.0 - alpha1))*(p_rghEqnComp2 & p_rgh)/rho_foam
              - pos(alpha1)*(p_rghEqnComp1 & p_rgh)/rho_gas
            );

            phi = phiHbyA + p_rghEqnIncomp.flux();
            U = HbyA
            + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }
    p = max(p_rgh + (alpha1*rho_gas + (1.0 - alpha1)*rho_foam)*gh, pMin);

    // Update densities from change in p_rgh
    rho_gas += Psi1*(p_rgh - p_rgh_0);
    rho_foam += Psi2*(p_rgh - p_rgh_0);
    rho = alpha1*rho_gas + (1 - alpha1)*rho_foam;
    K = 0.5*magSqr(U);
    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
