if (PBESwitch)
{
    if (blowingAgent == "no")
    {
        fvScalarMatrix wCO2_lEqn
        (
            fvm::ddt(wCO2_l)
          + fvm::div(phi,wCO2_l)
          + fvm::SuSp(expImpSource,wCO2_l)
          - fvm::div(alpha1phir,wCO2_l)
         == wCO2_lSource
        );
        wCO2_lEqn.relax();
        wCO2_lEqn.solve();

        fvScalarMatrix wCO2_gEqn
        (
            fvm::ddt(wCO2_g)
          + fvm::div(phi,wCO2_g)
          + fvm::SuSp(expImpSource,wCO2_g)
          - fvm::div(alpha1phir,wCO2_g)
         == wCO2_gSource
        );
        wCO2_gEqn.relax();
        wCO2_gEqn.solve();
    }
    else
    {
        fvScalarMatrix wBA_lEqn
        (
            fvm::ddt(wBA_l)
          + fvm::div(phi,wBA_l)
          + fvm::SuSp(expImpSource,wBA_l)
          - fvm::div(alpha1phir,wBA_l)
         == wBA_lSource
        );
        wBA_lEqn.relax();
        wBA_lEqn.solve();

        fvScalarMatrix wBA_gEqn
        (
            fvm::ddt(wBA_g)
          + fvm::div(phi,wBA_g)
          + fvm::SuSp(expImpSource,wBA_g)
          - fvm::div(alpha1phir,wBA_g)
         == wBA_gSource
        );
        wBA_gEqn.relax();
        wBA_gEqn.solve();

        fvScalarMatrix wCO2_lEqn
        (
            fvm::ddt(wCO2_l)
          + fvm::div(phi,wCO2_l)
          + fvm::SuSp(expImpSource,wCO2_l)
          - fvm::div(alpha1phir,wCO2_l)
         == wCO2_lSource
        );
        wCO2_lEqn.relax();
        wCO2_lEqn.solve();

        fvScalarMatrix wCO2_gEqn
        (
            fvm::ddt(wCO2_g)
          + fvm::div(phi,wCO2_g)
          + fvm::SuSp(expImpSource,wCO2_g)
          - fvm::div(alpha1phir,wCO2_g)
         == wCO2_gSource
        );
        wCO2_gEqn.relax();
        wCO2_gEqn.solve();
    }
}
else
{
    if (blowingAgent == "no")
    {
        fvScalarMatrix wCO2_gEqn
        (
            fvm::ddt(wCO2_g)
          + fvm::div(phi,wCO2_g)
          + fvm::SuSp(expImpSource,wCO2_g)
          - fvm::div(alpha1phir,wCO2_g)
         == wCO2_gSource
        );
        wCO2_gEqn.relax();
        wCO2_gEqn.solve();
    }
    else
    {
        fvScalarMatrix wBA_lEqn
        (
            fvm::ddt(wBA_l)
          + fvm::div(phi,wBA_l)
          + fvm::SuSp(expImpSource,wBA_l)
          - fvm::div(alpha1phir,wBA_l)
         == wBA_lSource
        );
        wBA_lEqn.relax();
        wBA_lEqn.solve();

        fvScalarMatrix wBA_gEqn
        (
            fvm::ddt(wBA_g)
          + fvm::div(phi,wBA_g)
          + fvm::SuSp(expImpSource,wBA_g)
          - fvm::div(alpha1phir,wBA_g)
         == wBA_gSource
        );
        wBA_gEqn.relax();
        wBA_gEqn.solve();

        fvScalarMatrix wCO2_gEqn
        (
            fvm::ddt(wCO2_g)
          + fvm::div(phi,wCO2_g)
          + fvm::SuSp(expImpSource,wCO2_g)
          - fvm::div(alpha1phir,wCO2_g)
         == wCO2_gSource
        );
        wCO2_gEqn.relax();
        wCO2_gEqn.solve();
    }
}