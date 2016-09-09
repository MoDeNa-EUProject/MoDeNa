{

    fvScalarMatrix TSEqn
    (
        fvm::ddt(TS)
      + fvm::div(phi,TS)
      + fvm::SuSp(expImpSource,TS)
      - fvm::div(alpha1phir,TS)
      - fvm::laplacian (thermalDiffusivity, TS)
     == TSSource
    );

    TSEqn.relax();
    TSEqn.solve();
}
