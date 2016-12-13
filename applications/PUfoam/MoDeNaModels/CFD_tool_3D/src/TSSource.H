double dXOHdt,XOH_C,XOH_Old,dWdt,XW_C,XW_Old,dLdt,alphaFO,
    rhoFO,C_TOT, L_l_C, L_l_Old, T_C,T_O, wBA_g_val, xBL_val, rhoPolyMixture;

// Calling the model for density reaction mixture
size_t Tm_pos =
    modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
size_t XOHv_pos =
    modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
modena_model_argPos_check(density_reaction_mixturemodel);

forAll(mesh.C(), celli)
{
    alphaFO = alpha2[celli];
    rhoFO = rho_foam[celli];

    XOH_C = XOH[celli];
    XOH_Old = XOH.oldTime()[celli];
    dXOHdt = (mag(XOH_C-XOH_Old))/(max(runTime.deltaTValue(),ROOTVSMALL));

    XW_C = XW[celli];
    XW_Old = XW.oldTime()[celli];
    dWdt = (mag(XW_C-XW_Old))/(max(runTime.deltaTValue(),ROOTVSMALL));

    T_C = TS[celli];
    T_O = TS.oldTime()[celli];

    wBA_g_val = wBA_g[celli];

    if (PBESwitch)
    {
        L_l_C = wBA_l[celli];
        L_l_Old = wBA_l.oldTime().oldTime()[celli];
        dLdt = (L_l_C - L_l_Old)/(max(runTime.deltaTValue(),ROOTVSMALL));
    }
    else
    {
        double lmax;
        lmax    = LliqMax(T_C);
        if (lmax < L0)
        {
            dLdt =
            (
                -1.0*ddT_LliqMax(T_C)*((T_C-T_O)
                /(max(runTime.deltaTValue(),ROOTVSMALL)))
            );
        }
        else
        {
            dLdt = 0.0;
        }
    }

    if (liquidMixtureDensitySurrogate)
    {
        // set input vector
        modena_inputs_set(inputs_den, Tm_pos, T_C);
        modena_inputs_set(inputs_den, XOHv_pos, XOH_C);
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
        rhoPolyMixture = modena_outputs_get(outputs_den, 0);
    }
    else
    {
        rhoPolyMixture = rhoPoly;
    }

    if (blowingAgent == "n-pentane")
    {
        C_TOT = C_Poly;
        if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
        {
            TSSource[celli] =
            (
                ((-DH_OH*COH_0*dXOHdt)/(rhoPolyMixture*C_TOT)
              + (-DH_W*CW_0*dWdt)/(rhoPolyMixture*C_TOT)
              + (latenth*dLdt)/(C_TOT))
            );
        }
        else if (alphaFO > 0.5 && XOH[celli] >= XOH_Gel)
        {
            TSSource[celli] = (latenth*dLdt)/(C_TOT);
        }
        else
        {
            TSSource[celli] = ROOTVSMALL;
        }
    }

    if (blowingAgent == "R-11")
    {
        if (wBA_g_val > 0.0)
        {
            xBL_val = xBL(T_C, dxdT);

            C_TOT =
            (
                C_Poly + (-(M_B/M_NCO)*
                (1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*latenth
            );

            if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
            {
                TSSource[celli] =
                (
                    ((-DH_OH*COH_0*dXOHdt)/(rhoPolyMixture*C_TOT)
                  + (-DH_W*CW_0*dWdt)/(rhoPolyMixture*C_TOT))
                );
            }
            else if (alphaFO > 0.5 && XOH[celli] >= XOH_Gel)
            {
                TSSource[celli] =
                (
                    ((-(M_B/M_NCO)*
                    (1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*latenth)
                );
            }
            else
            {
                TSSource[celli] = ROOTVSMALL;
            }
        }
        else
        {
            C_TOT = C_Poly;
            if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
            {
                TSSource[celli] =
                (
                    ((-DH_OH*COH_0*dXOHdt)/(rhoPolyMixture*C_TOT)
                  + (-DH_W*CW_0*dWdt)/(rhoPolyMixture*C_TOT))
                );
            }
            else
            {
                TSSource[celli] = ROOTVSMALL;
            }
        }
    }
    if (blowingAgent == "no")
    {
        C_TOT = C_Poly;
        if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
        {
            TSSource[celli] =
            (
                ((-DH_OH*COH_0*dXOHdt)/(rhoPolyMixture*C_TOT)
              + (-DH_W*CW_0*dWdt)/(rhoPolyMixture*C_TOT))
            );
        }
        else
        {
            TSSource[celli] = ROOTVSMALL;
        }
    }
    if (thermalConductivitySurrogateSwitch)
    {
        if (alphaFO > 0.5)
        {
            // surrogate model for thermal conductivity for CO2, cyclo-pentane
            modena_inputs_set(inputs_thermalConductivity, porosity_Pos, max(1.0 - rhoFO/rhoPolyMixture,0.0));
            scalar radiusBubble;
            radiusBubble = bubbleRadius(mZero[celli], mOne[celli]);
            modena_inputs_set(inputs_thermalConductivity, cell_size_Pos, (2.0*radiusBubble));
            modena_inputs_set(inputs_thermalConductivity, temp_Pos, TS[celli]);
            scalar pp_CO2, pp_BA;
            pp_CO2 =
                    partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension, wCO2_g[celli],
                        wBA_g[celli], p[celli], radiusBubble
                    );
            pp_BA =
               partialPressureBA
                        (
                            M_B, M_CO2, surfaceTension, wBA_g[celli], wCO2_g[celli],
                            p[celli], radiusBubble
                        );
            modena_inputs_set(inputs_thermalConductivity, X_CO2_Pos, (pp_CO2/(pp_BA+pp_CO2)));
            modena_inputs_set(inputs_thermalConductivity, X_O2_Pos, 0.0);
            modena_inputs_set(inputs_thermalConductivity, X_N2_Pos, 0.0);
            modena_inputs_set(inputs_thermalConductivity, X_Cyp_Pos, (pp_BA/(pp_BA+pp_CO2)));
            // input for strut content
            modena_inputs_set(inputs_strutContent, rho_foam_Pos, rhoFO);
            int ret_strutContent = modena_model_call (strutContentmodel, inputs_strutContent, outputs_strutContent);
            scalar st_c;
            st_c = modena_outputs_get(outputs_strutContent, 0);
            modena_inputs_set(inputs_thermalConductivity, strut_c_Pos, st_c);
            int ret_thermalConductivitymodel = modena_model_call (thermalConductivitymodel, inputs_thermalConductivity, outputs_thermalConductivity);
            if(modena_error_occurred())
            {
                exit(modena_error());
            }
            thermalConductivity[celli] = modena_outputs_get(outputs_thermalConductivity, 0);
            thermalDiff_foam[celli] = (thermalConductivity[celli]/(rhoFO*C_Poly));
        }
        else
        {
            thermalDiff_foam[celli] = thermalDiffusivityGas(T_C);
        }
        thermalDiff_gas[celli] = thermalDiffusivityGas(T_C);
    }
    else
    {
        if (blowingAgent == "n-pentane")
        {
            if (alphaFO > 0.5)
            {
                if (rhoFO > 48)
                {
                    thermalConductivity[celli] =
                    (
                        (scalar(8.7006e-8)*Foam::pow(rhoFO,2)
                      + scalar(8.4674e-5)*rhoFO + scalar(1.16e-2))
                    );
                }
                else
                {
                    thermalConductivity[celli] =
                    (
                        (scalar(9.3738e-6)*Foam::pow(rhoFO,2)
                      - scalar(7.3511e-4)*rhoFO+scalar(2.956e-2))
                    );
                }
                thermalDiff_foam[celli] = (thermalConductivity[celli]/(rhoFO*C_Poly));
            }
            else
            {
                thermalDiff_foam[celli] = thermalDiffusivityGas(T_C);
            }

            thermalDiff_gas[celli] = thermalDiffusivityGas(T_C);
        }
        else
        {
            if (alphaFO > 0.5)
            {
                thermalConductivity[celli]  = 0.03;
                thermalDiff_foam[celli] =
                    (thermalConductivity[celli]/(rhoFO*C_Poly));
            }
            else
            {
                thermalDiff_foam[celli] = thermalDiffusivityGas(T_C);
            }
            thermalDiff_gas[celli] = thermalDiffusivityGas(T_C);
        }
    }
    thermalDiffusivity[celli] =
        alpha1[celli]*thermalDiff_gas[celli]
      + alpha2[celli]*thermalDiff_foam[celli];
}
thermalDiffusivity.correctBoundaryConditions();
