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
if (KineticsModel == "simple")
{
    fvScalarMatrix EG_NCOEqn
    (
          fvm::ddt(rhoFoamAlpha2,EG_NCO)
        + fvm::div(frhoFoamAlpha2,EG_NCO)
      == source_EG_NCO
    );
    EG_NCOEqn.relax();
    EG_NCOEqn.solve();

    fvScalarMatrix EG_OHEqn
    (
        fvm::ddt(rhoFoamAlpha2,EG_OH)
        + fvm::div(frhoFoamAlpha2,EG_OH)
     == source_EG_OH
    );
    EG_OHEqn.relax();
    EG_OHEqn.solve();

    fvScalarMatrix H2OEqn
    (
        fvm::ddt(rhoFoamAlpha2,H2O)
      + fvm::div(frhoFoamAlpha2,H2O)
    == source_H2O
    );
    H2OEqn.relax();
    H2OEqn.solve();

    fvScalarMatrix CO2_SurEqn
    (
        fvm::ddt(rhoFoamAlpha2,CO2_Sur)
      + fvm::div(frhoFoamAlpha2,CO2_Sur)
    == source_CO2
    );
    CO2_SurEqn.relax();
    CO2_SurEqn.solve();

    fvScalarMatrix PENTANEEqn
    (
        fvm::ddt(rhoFoamAlpha2,PENTANE)
      + fvm::div(frhoFoamAlpha2,PENTANE)
     == source_PENTANE
    );
    PENTANEEqn.relax();
    PENTANEEqn.solve();

    fvScalarMatrix POLYMEREqn
    (
        fvm::ddt(rhoFoamAlpha2,POLYMER)
      + fvm::div(frhoFoamAlpha2,POLYMER)
     == source_POLYMER
    );
    POLYMEREqn.relax();
    POLYMEREqn.solve();

    fvScalarMatrix POLMERBLOWEqn
    (
        fvm::ddt(rhoFoamAlpha2,POLMERBLOW)
      + fvm::div(frhoFoamAlpha2,POLMERBLOW)
     == source_POLMERBLOW
    );
    POLMERBLOWEqn.relax();
    POLMERBLOWEqn.solve();

    fvScalarMatrix UREAEqn
    (
        fvm::ddt(rhoFoamAlpha2,UREA)
      + fvm::div(frhoFoamAlpha2,UREA)
     == source_UREA
    );
    UREAEqn.relax();
    UREAEqn.solve();

    fvScalarMatrix R_1_tempEqn
    (
        fvm::ddt(rhoFoamAlpha2,R_1_temp)
      + fvm::div(frhoFoamAlpha2,R_1_temp)
     == source_R_1_temp
    );
    R_1_tempEqn.relax();
    R_1_tempEqn.solve();
}
if (KineticsModel == "RF-1")
{
    fvScalarMatrix Catalyst_1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,Catalyst_1)
      + fvm::div(frhoFoamAlpha2,Catalyst_1)
     == source_Catalyst_1
    );
    Catalyst_1Eqn.relax();
    Catalyst_1Eqn.solve();

    fvScalarMatrix CE_A0Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_A0)
      + fvm::div(frhoFoamAlpha2,CE_A0)
     == source_CE_A0 
    );
    CE_A0Eqn.relax();
    CE_A0Eqn.solve();

    fvScalarMatrix CE_A1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_A1)
      + fvm::div(frhoFoamAlpha2,CE_A1)
     == source_CE_A1
    );
    CE_A1Eqn.relax();
    CE_A1Eqn.solve();

    fvScalarMatrix CE_BEqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_B)
      + fvm::div(frhoFoamAlpha2,CE_B)
     == source_CE_B
    );
    CE_BEqn.relax();
    CE_BEqn.solve();

    fvScalarMatrix CE_B2Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_B2)
      + fvm::div(frhoFoamAlpha2,CE_B2)
     == source_CE_B2
    );
    CE_B2Eqn.relax();
    CE_B2Eqn.solve();

    fvScalarMatrix CE_I0Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_I0)
      + fvm::div(frhoFoamAlpha2,CE_I0)
     == source_CE_I0
    );
    CE_I0Eqn.relax();
    CE_I0Eqn.solve();

    fvScalarMatrix CE_I1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_I1)
      + fvm::div(frhoFoamAlpha2,CE_I1)
     == source_CE_I1
    );
    CE_I1Eqn.relax();
    CE_I1Eqn.solve();

    fvScalarMatrix CE_I2Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_I2)
      + fvm::div(frhoFoamAlpha2,CE_I2)
     == source_CE_I2
    );
    CE_I2Eqn.relax();
    CE_I2Eqn.solve();

    fvScalarMatrix CE_PBAEqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_PBA)
      + fvm::div(frhoFoamAlpha2,CE_PBA)
     == source_CE_PBA
    );
    CE_PBAEqn.relax();
    CE_PBAEqn.solve();

    fvScalarMatrix CE_BreacEqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_Breac)
      + fvm::div(frhoFoamAlpha2,CE_Breac)
     == source_CE_Breac
    );
    CE_BreacEqn.relax();
    CE_BreacEqn.solve();

    fvScalarMatrix CE_Areac0Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_Areac0)
      + fvm::div(frhoFoamAlpha2,CE_Areac0)
     == source_CE_Areac0
    );
    CE_Areac0Eqn.relax();
    CE_Areac0Eqn.solve();

    fvScalarMatrix CE_Areac1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_Areac1)
      + fvm::div(frhoFoamAlpha2,CE_Areac1)
     == source_CE_Areac1
    );
    CE_Areac1Eqn.relax();
    CE_Areac1Eqn.solve();

    fvScalarMatrix CE_Ireac0Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_Ireac0)
      + fvm::div(frhoFoamAlpha2,CE_Ireac0)
     == source_CE_Ireac0
    );
    CE_Ireac0Eqn.relax();
    CE_Ireac0Eqn.solve();

    fvScalarMatrix CE_Ireac1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_Ireac1)
      + fvm::div(frhoFoamAlpha2,CE_Ireac1)
     == source_CE_Ireac1
    );
    CE_Ireac1Eqn.relax();
    CE_Ireac1Eqn.solve();

    fvScalarMatrix CE_Ireac2Eqn
    (
        fvm::ddt(rhoFoamAlpha2,CE_Ireac2)
      + fvm::div(frhoFoamAlpha2,CE_Ireac2)
     == source_CE_Ireac2
    );
    CE_Ireac2Eqn.relax();
    CE_Ireac2Eqn.solve();

    fvScalarMatrix BulkEqn
    (
        fvm::ddt(rhoFoamAlpha2,Bulk)
      + fvm::div(frhoFoamAlpha2,Bulk)
     == source_Bulk
    );
    BulkEqn.relax();
    BulkEqn.solve();

    fvScalarMatrix R_1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,R_1)
      + fvm::div(frhoFoamAlpha2,R_1)
     == source_R_1
    );
    R_1Eqn.relax();
    R_1Eqn.solve();

    fvScalarMatrix R_1_temp_RF1Eqn
    (
        fvm::ddt(rhoFoamAlpha2,R_1_temp_RF1)
      + fvm::div(frhoFoamAlpha2,R_1_temp_RF1)
     == source_R_1_temp_RF1
    );
    R_1_temp_RF1Eqn.relax();
    R_1_temp_RF1Eqn.solve();
}
