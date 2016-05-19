alpha1phir = fvc::interpolate(alpha1)*phir;

rhoFoamAlpha2  = rho_foam*alpha2;
frhoFoamAlpha2 = fvc::interpolate(rhoFoamAlpha2)*phi;

volScalarField expImpSource
(
    fvc::div(alpha1phir) - fvc::div(phi)
);

if (KineticsModel == "generic")
{

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
}

if (KineticsModel == "RF-1")
{

    fvScalarMatrix XWEqn
    (
        fvm::ddt(XW)
      + fvm::div(phi,XW)
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
     == XOHSource
    );
    XOHEqn.relax();
    XOHEqn.solve();

    fvScalarMatrix Catalyst_1Eqn
    (
        fvm::ddt(Catalyst_1)
      + fvm::div(phi,Catalyst_1)
      + fvm::SuSp(expImpSource,Catalyst_1)
      - fvm::div(alpha1phir,Catalyst_1)
     == source_Catalyst_1
    );
    Catalyst_1Eqn.relax();
    Catalyst_1Eqn.solve();

    fvScalarMatrix CE_A0Eqn
    (
        fvm::ddt(CE_A0)
      + fvm::div(phi,CE_A0)
      + fvm::SuSp(expImpSource,CE_A0)
      - fvm::div(alpha1phir,CE_A0)
     == source_CE_A0
    );
    CE_A0Eqn.relax();
    CE_A0Eqn.solve();

    fvScalarMatrix CE_A1Eqn
    (
        fvm::ddt(CE_A1)
      + fvm::div(phi,CE_A1)
      + fvm::SuSp(expImpSource,CE_A1)
      - fvm::div(alpha1phir,CE_A1)
     == source_CE_A1
    );
    CE_A1Eqn.relax();
    CE_A1Eqn.solve();

    fvScalarMatrix CE_BEqn
    (
        fvm::ddt(CE_B)
      + fvm::div(phi,CE_B)
      + fvm::SuSp(expImpSource,CE_B)
      - fvm::div(alpha1phir,CE_B)
     == source_CE_B
    );
    CE_BEqn.relax();
    CE_BEqn.solve();

    fvScalarMatrix CE_B2Eqn
    (
        fvm::ddt(CE_B2)
      + fvm::div(phi,CE_B2)
      + fvm::SuSp(expImpSource,CE_B2)
      - fvm::div(alpha1phir,CE_B2)
     == source_CE_B2
    );
    CE_B2Eqn.relax();
    CE_B2Eqn.solve();

    fvScalarMatrix CE_I0Eqn
    (
        fvm::ddt(CE_I0)
      + fvm::div(phi,CE_I0)
      + fvm::SuSp(expImpSource,CE_I0)
      - fvm::div(alpha1phir,CE_I0)
      == source_CE_I0
    );
    CE_I0Eqn.relax();
    CE_I0Eqn.solve();

    fvScalarMatrix CE_I1Eqn
    (
        fvm::ddt(CE_I1)
      + fvm::div(phi,CE_I1)
      + fvm::SuSp(expImpSource,CE_I1)
      - fvm::div(alpha1phir,CE_I1)
      == source_CE_I1
    );
    CE_I1Eqn.relax();
    CE_I1Eqn.solve();

    fvScalarMatrix CE_I2Eqn
    (
        fvm::ddt(CE_I2)
      + fvm::div(phi,CE_I2)
      + fvm::SuSp(expImpSource,CE_I2)
      - fvm::div(alpha1phir,CE_I2)
     == source_CE_I2
    );
    CE_I2Eqn.relax();
    CE_I2Eqn.solve();

    fvScalarMatrix CE_PBAEqn
    (
        fvm::ddt(CE_PBA)
      + fvm::div(phi,CE_PBA)
      + fvm::SuSp(expImpSource,CE_PBA)
      - fvm::div(alpha1phir,CE_PBA)
     == source_CE_PBA
    );
    CE_PBAEqn.relax();
    CE_PBAEqn.solve();

    fvScalarMatrix CE_BreacEqn
    (
        fvm::ddt(CE_Breac)
      + fvm::div(phi,CE_Breac)
      + fvm::SuSp(expImpSource,CE_Breac)
      - fvm::div(alpha1phir,CE_Breac)
     == source_CE_Breac
    );
    CE_BreacEqn.relax();
    CE_BreacEqn.solve();

    fvScalarMatrix CE_Areac0Eqn
    (
        fvm::ddt(CE_Areac0)
      + fvm::div(phi,CE_Areac0)
      + fvm::SuSp(expImpSource,CE_Areac0)
      - fvm::div(alpha1phir,CE_Areac0)
     == source_CE_Areac0
    );
    CE_Areac0Eqn.relax();
    CE_Areac0Eqn.solve();

    fvScalarMatrix CE_Areac1Eqn
    (
        fvm::ddt(CE_Areac1)
      + fvm::div(phi,CE_Areac1)
      + fvm::SuSp(expImpSource,CE_Areac1)
      - fvm::div(alpha1phir,CE_Areac1)
     == source_CE_Areac1
    );
    CE_Areac1Eqn.relax();
    CE_Areac1Eqn.solve();

    fvScalarMatrix CE_Ireac0Eqn
    (
        fvm::ddt(CE_Ireac0)
      + fvm::div(phi,CE_Ireac0)
      + fvm::SuSp(expImpSource,CE_Ireac0)
      - fvm::div(alpha1phir,CE_Ireac0)
     == source_CE_Ireac0
    );
    CE_Ireac0Eqn.relax();
    CE_Ireac0Eqn.solve();

    fvScalarMatrix CE_Ireac1Eqn
    (
        fvm::ddt(CE_Ireac1)
      + fvm::div(phi,CE_Ireac1)
      + fvm::SuSp(expImpSource,CE_Ireac1)
      - fvm::div(alpha1phir,CE_Ireac1)
     == source_CE_Ireac1
    );
    CE_Ireac1Eqn.relax();
    CE_Ireac1Eqn.solve();

    fvScalarMatrix CE_Ireac2Eqn
    (
        fvm::ddt(CE_Ireac2)
      + fvm::div(phi,CE_Ireac2)
      + fvm::SuSp(expImpSource,CE_Ireac2)
      - fvm::div(alpha1phir,CE_Ireac2)
     == source_CE_Ireac2
    );
    CE_Ireac2Eqn.relax();
    CE_Ireac2Eqn.solve();

    fvScalarMatrix BulkEqn
    (
        fvm::ddt(Bulk)
      + fvm::div(phi,Bulk)
      + fvm::SuSp(expImpSource,Bulk)
      - fvm::div(alpha1phir,Bulk)
     == source_Bulk
    );
    BulkEqn.relax();
    BulkEqn.solve();

    fvScalarMatrix R_1Eqn
    (
        fvm::ddt(R_1)
      + fvm::div(phi,R_1)
      + fvm::SuSp(expImpSource,R_1)
      - fvm::div(alpha1phir,R_1)
     == source_R_1
    );
    R_1Eqn.relax();
    R_1Eqn.solve();

    fvScalarMatrix R_1_temp_RF1Eqn
    (
        fvm::ddt(R_1_temp_RF1)
      + fvm::div(phi,R_1_temp_RF1)
      + fvm::SuSp(expImpSource,R_1_temp_RF1)
      - fvm::div(alpha1phir,R_1_temp_RF1)
     == source_R_1_temp_RF1
    );
    R_1_temp_RF1Eqn.relax();
    R_1_temp_RF1Eqn.solve();
}
