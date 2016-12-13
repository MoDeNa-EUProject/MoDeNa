if (PBESwitch)
{
    scalar prsr, tmpt, alpha_f, g1_BA_val,
           cc1_value, tmptrOld, timeStep, dTdt,
           dXWdt, g1_CO2_val, XOH_val, rhoPolyS;
    scalar henry_Coefficient, mZero_val, mOne_val,
           bubble_radius, partialPressure_CO2, wCO2_g_val,
           wBA_g_val, wCO2_Max,wBL_Max, partialPressure_BA;

    // call surrogate model for density reaction mixture
    size_t T_pos =
        modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
    size_t XOH_pos =
        modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
    modena_model_argPos_check(density_reaction_mixturemodel);

    forAll(mesh.C(), celli)
    {
        tmpt = TS[celli];
        tmptrOld = TS.oldTime()[celli];
        // tmptrOldOld = TS.oldTime().oldTime()[celli];
        timeStep = runTime.deltaTValue();
        prsr = p[celli];
        // rho_f = rho_foam[celli];
        alpha_f = alpha2[celli];
        mZero_val = mZero[celli];
        mOne_val = mOne[celli];
        wCO2_g_val = wCO2_g[celli];
        wBA_g_val = wBA_g[celli];
        cc1_value = cc1[celli];
        g1_BA_val = g1_BA[celli];
        XOH_val = XOH[celli];

        if (liquidMixtureDensitySurrogate)
        {
            // set input vector
            modena_inputs_set(inputs_den, T_pos, tmpt);
            modena_inputs_set(inputs_den, XOH_pos, XOH_val);
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

        bubble_radius = bubbleRadius(mZero_val, mOne_val);

        partialPressure_BA = partialPressureBA
                        (
                            M_B, M_CO2, surfaceTension, wBA_g_val,
                            wCO2_g_val, prsr, bubble_radius
                        );
        if (blowingAgent == "n-pentane")
        {
            if (alpha_f > 0.5)
            {
                // if (LliqMax(tmpt) < L0)
                // {
                    // dTdt  = dTdtFirstOrder(tmpt, tmptrOld, timeStep);
                    // wBA_gSource[celli] = (-cc1_value*ddT_LliqMax(tmpt)*dTdt);

                    wBA_lSource[celli] =
                        (
                            (M_B/1000.0)*(1.0/rhoPolyS)
                            *(partialPressure_BA/(max((RR*tmpt),ROOTVSMALL)))*(-g1_BA_val)
                        );
                        wBA_gSource[celli] = -wBA_lSource[celli];
                // }
                // else
                // {
                    // wBA_gSource[celli] = 0.0;
                    // wBA_lSource[celli] = 0.0;
                // }
            }
            else
            {
                wBA_lSource[celli] = 0.0;
                wBA_gSource[celli] = 0.0;
            }
        }

        if (blowingAgent == "R-11")
        {
            scalar xBL_val;
            xBL_val = xBL(tmpt, dxdT);
            wBL_Max = wBL_D(xBL_val, M_B, M_NCO, L0);
            if (wBL_Max > 1.0e-4 && alpha_f > 0.5)
            {
                dTdt = dTdtFirstOrder(tmpt, tmptrOld, timeStep);
                wBA_gSource[celli] =
                    (
                        cc1_value*(-(M_B/M_NCO)
                        *(1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*dTdt
                    );
                wBA_lSource[celli] =
                    (
                        -cc1_value*(-(M_B/M_NCO)
                        *(1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*dTdt
                    );
            }
            else
            {
                wBA_gSource[celli] = 0.0;
                wBA_lSource[celli] = 0.0;
            }
        }
        if (alpha_f > 0.5)
        {
            g1_CO2_val = g1_CO2[celli];
            dXWdt =
                (
                    (mag(XW[celli] - XW.oldTime()[celli]))
                    /(max(runTime.deltaTValue(),ROOTVSMALL))
                );
            partialPressure_CO2 =
                partialPressureCO2
                (
                    M_CO2, M_B, surfaceTension,
                    wCO2_g_val, wBA_g_val, prsr, bubble_radius
                );
            if (CW_0 == 0.0)
            {
                wCO2_lSource[celli] = 0.0;
            }
            else
            {
                if (XOH[celli] < XOH_Gel)
                {
                    wCO2_lSource[celli] =
                    (
                        (CW_0*dXWdt*(M_CO2/1000.0)*(1.0/rhoPolyS)
                        - g1_CO2_val*(partialPressure_CO2/(RR*tmpt))*(M_CO2/1000.0)*(1.0/rhoPolyS))
                    );
                }
            }

            henry_Coefficient = henryCoefficient(tmpt);

            wCO2_Max = wCO2Max(M_CO2, M_liq, partialPressure_CO2, henry_Coefficient);

            if (wCO2_l[celli] < wCO2_Max || CW_0 == 0.0 || XOH[celli] > XOH_Gel)
            {
                wCO2_gSource[celli] = 0.0;

            }
            else
            {
                wCO2_gSource[celli] =
                    (
                        g1_CO2_val
                        *(partialPressure_CO2/(RR*tmpt))*(M_CO2/1000.0)*(1.0/rhoPolyS)
                    );
            }
        }
        else
        {
            wCO2_gSource[celli] = 0.0;
            wCO2_lSource[celli] = 0.0;
        }
    }
}
else
{
    scalar prsr, tmpt, alpha_f,
           QBA, tmptrOld, timeStep, dTdt;

    scalar henry_Coefficient,
           bubble_radius, partialPressure_CO2, wCO2_g_val,
           wBA_g_val, wCO2_Max, wBL_Max, XOH_val, rhoPolyS;

    scalar m0_val, m1_val;
    m0_val = m1_val = 0.0;
    scalar xBL_val;

    // call surrogate model for density reaction mixture
    size_t T_pos =
        modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
    size_t XOH_pos =
        modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
    modena_model_argPos_check(density_reaction_mixturemodel);

    forAll(mesh.C(), celli)
    {
        tmpt = TS[celli];
        tmptrOld = TS.oldTime()[celli];
        // tmptOldOld = TS.oldTime().oldTime()[celli];
        timeStep = runTime.deltaTValue();
        prsr = p[celli];
        // rho_f = rho_foam[celli];
        alpha_f = alpha2[celli];
        wCO2_g_val = wCO2_g[celli];
        wBA_g_val = wBA_g[celli];
        // wBA_l_val = wBA_l[celli];
        XOH_val = XOH[celli];

        if (liquidMixtureDensitySurrogate)
        {
            // set input vector
            modena_inputs_set(inputs_den, T_pos, tmpt);
            modena_inputs_set(inputs_den, XOH_pos, XOH_val);
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

        if (blowingAgent == "n-pentane")
        {
            if (alpha_f > 0.5)
            {
                dTdt = dTdtFirstOrder(tmpt, tmptrOld, timeStep);

                if (L0 < LliqMax(tmpt))
                {
                    QBA = -ddT_LliqMax(tmpt)*dTdt;
                    wBA_gSource[celli] = QBA;
                    wBA_lSource[celli] = ddT_LliqMax(tmpt)*dTdt;
                }
                else
                {
                    wBA_gSource[celli] = 0.0;
                    wBA_lSource[celli] = 0.0;
                }
            }
            else
            {
                wBA_lSource[celli] = 0.0;
                wBA_gSource[celli] = 0.0;
            }
        }

        if (blowingAgent == "R-11")
        {
            xBL_val = xBL(tmpt, dxdT);
            wBL_Max = wBL_D(xBL_val, M_B, M_NCO, L0);
            if (wBL_Max > 1.0e-4 && alpha_f > 0.5)
            {
                dTdt = dTdtFirstOrder(tmpt, tmptrOld, timeStep);
                wBA_gSource[celli] =
                    (
                        (-(M_B/M_NCO)
                        *(1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*dTdt
                    );
                wBA_lSource[celli] =
                    (
                        ((M_B/M_NCO)
                        *(1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*dTdt
                    );
            }
            else
            {
                wBA_lSource[celli] = 0.0;
                wBA_gSource[celli] = 0.0;
            }
        }

        wCO2_lSource[celli] = scalar(0.0);

        henry_Coefficient = henryCoefficient(tmpt);
        bubble_radius = bubbleRadius(m0_val, m1_val);
        partialPressure_CO2 =
            partialPressureCO2
            (
                M_CO2, M_B, surfaceTension,
                wCO2_g_val, wBA_g_val, prsr, bubble_radius
            );
        wCO2_Max = wCO2Max(M_CO2, M_liq, partialPressure_CO2, henry_Coefficient);

        if (alpha_f > 0.5 && XOH[celli] <= XOH_Gel && CW_0 != 0.0)
        {
            wCO2_gSource[celli] =
            (
                ((CW_0*XW[celli]*(M_CO2/1000.0)*(1.0/rhoPolyS)) - wCO2_Max)
            );
        }
        else
        {
            wCO2_gSource[celli] = 0.0;
        }

        if (wCO2_gSource[celli] < 0.0)
        {
            wCO2_gSource[celli] = scalar(0.0);
        }
    }
}
