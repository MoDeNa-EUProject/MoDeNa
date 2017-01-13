/**
@ingroup mod_3Dcfd
@file TSEqn.H
@brief builds and solves the temperature equation
*/
///@cond
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
///@endcond