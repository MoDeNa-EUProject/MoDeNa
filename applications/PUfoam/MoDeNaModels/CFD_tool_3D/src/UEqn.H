volTensorField gradU = fvc::grad(U);
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U) 
  - fvm::laplacian(muMixture, U) 
  // - (fvc::grad(U) & fvc::grad(muMixture))
);

UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    K = 0.5*magSqr(U);
}