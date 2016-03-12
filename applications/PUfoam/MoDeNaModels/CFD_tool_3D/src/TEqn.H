{

volScalarField rhoAlphabar=rho*thermalDiffusivity;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T) 
      - fvm::laplacian (rhoAlphabar, T)
    == TSource
    );

    TEqn.relax();
    TEqn.solve();

    mixture.correct();
}
