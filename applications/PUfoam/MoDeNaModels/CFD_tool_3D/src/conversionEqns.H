alpha1phir = fvc::interpolate(alpha1)*phir;

volScalarField expImpSource
(
    fvc::div(alpha1phir) - fvc::div(phi)
);

fvScalarMatrix XWEqn
(
    fvm::ddt(XW)
  + fvm::div(phi,XW)
  + fvm::SuSp(XWSource,XW)
  + fvm::SuSp(expImpSource,XW)
  - fvm::div(alpha1phir,XW)
 == XWSource    
);
XWEqn.relax();
XWEqn.solve();

fvScalarMatrix XOHEqn
(
    fvm::ddt(XOH)
  + fvm::div(phi,XOH)
  + fvm::SuSp(expImpSource,XOH)
  - fvm::div(alpha1phir,XOH)
==
  - fvm::SuSp(Coeff_1stOrd,XOH)
  + ScNonXOH + SuXOH + fvm::Sp(SpXOH,XOH)
);
XOHEqn.relax();
XOHEqn.solve();

