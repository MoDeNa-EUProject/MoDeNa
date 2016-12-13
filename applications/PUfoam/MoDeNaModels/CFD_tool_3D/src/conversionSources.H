if (KineticsModel == "generic")
{
    scalar temp_current,
           temp_old,
           Arhns_current,
           Arhns_old,
           XW_val,
           XOH_valOld,
           BA_dilution,
           XOHSource_Old,
           QW,
           rhoPolyS;

    // call surrogate model for density reaction mixture
        size_t T_pos =
            modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
        size_t XOH_pos =
            modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
        modena_model_argPos_check(density_reaction_mixturemodel);

    forAll(mesh.C(), celli)
    {
        temp_current = TS[celli];
        temp_old = TS.oldTime()[celli];

        if (liquidMixtureDensitySurrogate)
        {
            // set input vector
            modena_inputs_set(inputs_den, T_pos, temp_current);
            modena_inputs_set(inputs_den, XOH_pos, XOH[celli]);
            // call the model
            int ret_den =
                modena_model_call
                (
                    density_reaction_mixturemodel,
                    inputs_den,
                    outputs_den
                );
            if (ret_den != 0)
                {
                    modena_inputs_destroy (inputs_den);
                    modena_outputs_destroy (outputs_den);
                    modena_model_destroy (density_reaction_mixturemodel);
                    exit(ret_den);
                }
            rhoPolyS = modena_outputs_get(outputs_den, 0);
        }
        else
        {
            rhoPolyS = rhoPoly;
        }

        if (XW[celli] < 0.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        if (XW[celli] > 1.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        QW = QKinW (AW, EW, temp_current, L0, rhoPolyS, rhoBL);
        if (QW < 0.0 || CW_0 == 0.0)
        {
            QW = ROOTVSMALL;
        }

        if (XOH[celli] < XOH_Gel)
        {
            XWSource[celli] = QW;
        }

        if (XW[celli] < 0.0)
        {
            XW_val = ROOTVSMALL;
        }
        else if (XW[celli] > 1.0)
        {
            XW_val = ROOTVSMALL;
        }
        else
        {
            XW_val = XW[celli];
        }

        if (XOH[celli] < 0.0)
        {
            XOH[celli] = ROOTVSMALL;
        }
        if (XOH[celli] > 1.0)
        {
            XOH[celli] = ROOTVSMALL;
        }

        if (XOH[celli] < XOH_Gel)
        {
            XOH_valOld = XOH[celli];
            Arhns_current = arrhenius (AOH, EOH, temp_current);
            Arhns_old = arrhenius (AOH, EOH, temp_old);
            BA_dilution = BAdilution(L0,rhoPolyS,rhoBL);
            XOHSource_Old =
                (
                    (Arhns_old*COH_0*(1 - XOH_valOld)*
                    (CNCO_0/COH_0 - 2*XW_val*CW_0/COH_0 - XOH_valOld)*BA_dilution)
                );
            if (XOHSource_Old < 0.0)
            {
                ScNonXOH[celli] = ROOTVSMALL;
                Coeff_1stOrd[celli] = ROOTVSMALL;
                SuXOH[celli] = ROOTVSMALL;
                SpXOH[celli] = ROOTVSMALL;
            }
            else
            {
                ScNonXOH[celli] =
                    (
                        BA_dilution*
                        (Arhns_current*CNCO_0 - 2.0*Arhns_current*CW_0*XW_val)
                    );
                Coeff_1stOrd[celli] =
                (
                    BA_dilution*
                    (Arhns_current*COH_0 + Arhns_current*CNCO_0
                    - 2.0*Arhns_current*CW_0*XW_val)
                );

                SuXOH[celli] =
                (
                    -1.0*BA_dilution*Arhns_old*
                    COH_0*Foam::pow(XOH_valOld,2.0)
                );
                SpXOH[celli] =
                (
                    2.0*BA_dilution*Arhns_old*COH_0*XOH_valOld
                );
            }
        }
    }
}
if (KineticsModel == "RF-1")
{
    forAll(mesh.C(), celli)
    {
        modena_inputs_set(inputs_kinetics, kineticTime_Pos, runTime.value());
        modena_inputs_set(inputs_kinetics, Catalyst_1_Pos, (Catalyst_1[celli]));
        modena_inputs_set(inputs_kinetics, CE_A0_Pos, (CE_A0[celli]));
        modena_inputs_set(inputs_kinetics, CE_A1_Pos, (CE_A1[celli]));
        modena_inputs_set(inputs_kinetics, CE_B_Pos, (CE_B[celli]));
        modena_inputs_set(inputs_kinetics, CE_B2_Pos, (CE_B2[celli]));
        modena_inputs_set(inputs_kinetics, CE_I0_Pos, (CE_I0[celli]));
        modena_inputs_set(inputs_kinetics, CE_I1_Pos, (CE_I1[celli]));
        modena_inputs_set(inputs_kinetics, CE_I2_Pos, (CE_I2[celli]));
        modena_inputs_set(inputs_kinetics, CE_PBA_Pos, (CE_PBA[celli]));
        modena_inputs_set(inputs_kinetics, CE_Breac_Pos, (CE_Breac[celli]));
        modena_inputs_set(inputs_kinetics, CE_Areac0_Pos, (CE_Areac0[celli]));
        modena_inputs_set(inputs_kinetics, CE_Areac1_Pos, (CE_Areac1[celli]));
        modena_inputs_set(inputs_kinetics, CE_Ireac0_Pos, (CE_Ireac0[celli]));
        modena_inputs_set(inputs_kinetics, CE_Ireac1_Pos, (CE_Ireac1[celli]));
        modena_inputs_set(inputs_kinetics, CE_Ireac2_Pos, (CE_Ireac2[celli]));
        modena_inputs_set(inputs_kinetics, Bulk_Pos, Bulk[celli]);
        modena_inputs_set(inputs_kinetics, R_1_Pos, R_1[celli]);
        modena_inputs_set(inputs_kinetics, R_1_mass_Pos, R_1_mass[celli]);
        modena_inputs_set(inputs_kinetics, R_1_temp_RF1_Pos, (TS[celli] - 273.15));
        modena_inputs_set(inputs_kinetics, R_1_vol_RF1_Pos, R_1_vol_RF1[celli]);

        // call the model
        int ret_kinetics =
            modena_model_call(kinetics, inputs_kinetics, outputs_kinetics);
        if (ret_kinetics != 0)
        {
            modena_inputs_destroy (inputs_kinetics);
            modena_outputs_destroy (outputs_kinetics);
            modena_model_destroy (kinetics);
            exit(ret_kinetics);
        }
        // source terms for RF-1
        source_Catalyst_1[celli] =
            modena_outputs_get(outputs_kinetics, source_Catalyst_1_Pos);
        source_CE_A0[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_A0_Pos);
        source_CE_A1[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_A1_Pos);
        source_CE_B[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_B_Pos);
        source_CE_B2[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_B2_Pos);
        source_CE_I0[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_I0_Pos);
        source_CE_I1[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_I1_Pos);
        source_CE_I2[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_I2_Pos);
        source_CE_PBA[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_PBA_Pos);
        source_CE_Breac[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_Breac_Pos);
        source_CE_Areac0[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_Areac0_Pos);
        source_CE_Areac1[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_Areac1_Pos);
        source_CE_Ireac0[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_Ireac0_Pos);
        source_CE_Ireac1[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_Ireac1_Pos);
        source_CE_Ireac2[celli] =
            modena_outputs_get(outputs_kinetics, source_CE_Ireac2_Pos);
        source_Bulk[celli] =
            modena_outputs_get(outputs_kinetics, source_Bulk_Pos);
        source_R_1[celli] =
            modena_outputs_get(outputs_kinetics, source_R_1_Pos);
        source_R_1_mass[celli] =
            modena_outputs_get(outputs_kinetics, source_R_1_mass_Pos);
        source_R_1_temp_RF1[celli] =
            modena_outputs_get(outputs_kinetics, source_R_1_temp_RF1_Pos);
        source_R_1_vol_RF1[celli] =
            modena_outputs_get(outputs_kinetics, source_R_1_vol_RF1_Pos);

            if (CW_0 > scalar(1e-8)) {
                XWSource[celli] = scalar(-1.0)*
                    source_CE_B2[celli]/(CW_0/scalar(1000.0));
            } else {
                XWSource[celli] = scalar(0.0);
            }

            XOHSource[celli] = scalar(-1.0)*
                               (source_CE_A0[celli] + source_CE_A1[celli])/
                               ((COH_0))*scalar(1000.0);
    }
    source_Catalyst_1.correctBoundaryConditions();
    source_CE_A0.correctBoundaryConditions();
    source_CE_A1.correctBoundaryConditions();
    source_CE_B.correctBoundaryConditions();
    source_CE_B2.correctBoundaryConditions();
    source_CE_I0.correctBoundaryConditions();
    source_CE_I1.correctBoundaryConditions();
    source_CE_I2.correctBoundaryConditions();
    source_CE_PBA.correctBoundaryConditions();
    source_CE_Breac.correctBoundaryConditions();
    source_CE_Areac0.correctBoundaryConditions();
    source_CE_Areac1.correctBoundaryConditions();
    source_CE_Ireac0.correctBoundaryConditions();
    source_CE_Ireac1.correctBoundaryConditions();
    source_CE_Ireac2.correctBoundaryConditions();
    source_Bulk.correctBoundaryConditions();
    source_R_1.correctBoundaryConditions();
    source_R_1_mass.correctBoundaryConditions();
    source_R_1_temp_RF1.correctBoundaryConditions();
    source_R_1_vol_RF1.correctBoundaryConditions();
}
