if (PBESwitch)
{   
    if 
    (
        PBEMethod == "QMOM" 
     && nNodes == 2
    )
    {
        static int nodes = nNodes;
        int mOrder[2*nodes]; 
        for (int i = 0; i < 2*nNodes; i++){mOrder[i] = i;}
        double mom[2*nodes], we[nodes], vi[nodes], 
        sgBA[2*nodes], sgCO2[2*nodes], sc[2*nodes];
        double sources[2*nodes];

        double CO2_l_val, L_l_val, tmp_val, alphafoam;
        double cc1_val = 0.0;
        double LMAX, xBL_val, rhoPolyS;
        // Calling the model for density reaction mixture
        size_t T_pos = 
            modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
        size_t XOH_pos = 
            modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
        modena_model_argPos_check(density_reaction_mixturemodel);
        
        forAll(mesh.C(), celli)
        {
            CO2_l_val = wCO2_l[celli];
            L_l_val = wBA_l[celli];
            tmp_val = TS[celli]; 
            alphafoam = alpha2[celli];
            double wCO2_g_val, wBA_g_val,prsr;
            wCO2_g_val= wCO2_g[celli]; 
            wBA_g_val = wBA_g[celli];
            prsr = p[celli];

            if (liquidMixtureDensitySurrogate)
            {
                // set input vector
                modena_inputs_set(inputs_den, Td_pos, TS[celli]);
                modena_inputs_set(inputs_den, XOHd_pos, XOH[celli]);
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

            // Using old moments to compute the new weights and nodes
            if (alphafoam > 0.5)
            {
                mom[0] = mZero[celli];
                mom[1] = mOne[celli]; 
                mom[2] = mTwo[celli]; 
                mom[3] = mThree[celli]; 

                // calling the functions, PDA, growthSource, coalescenceSource
                PDA(we, vi, mom, nodes);    
                weight0[celli] = we[0];
                weight1[celli] = we[1];
                node0[celli] = vi[0];
                node1[celli] = vi[1]; 

                // Required variables for wCO2_Max
                double henry_Coefficient, bubble_radius, 
                    partialPressure_CO2, wCO2_Max;
                henry_Coefficient = henryCoefficient(tmp_val);
                bubble_radius = bubbleRadius(mom[0], mom[1]);
                partialPressure_CO2 = 
                    partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension, wCO2_g_val,
                        wBA_g_val, prsr, bubble_radius
                    );
                wCO2_Max = 
                    wCO2Max
                    (
                        M_CO2, M_liq, partialPressure_CO2, henry_Coefficient
                    );
            
                if (blowingAgent == "n-pentane")
                {
                    LMAX  = LliqMax (tmp_val);
                }
                if (blowingAgent == "R-11")
                {
                    xBL_val = xBL(tmp_val, dxdT);
                    LMAX = wBL_D(xBL_val, M_B, M_NCO, L0);
                }
                if (blowingAgent == "no")
                {
                    LMAX = 0.0;
                }
                growthSource
                    (
                        sgBA, sgCO2, we, vi, nodes, mOrder, CO2_l_val, L_l_val,
                        tmp_val, wCO2_Max, cc1_val, LMAX
                    );
            
                g1_CO2[celli] = sgCO2[1];
                g1_BA[celli] = sgBA[1];
                cc1[celli] = cc1_val;
            
                if (XOH[celli] > XOH_Gel)
                {
                    for (int j = 0; j < 2*nodes; j++)
                    {
                        sc[j] = 0.0;
                    }
                }   
                else
                {
                    coalescenceSource(sc, we, vi, nodes, mOrder);
                }
            
                for (int j = 0; j < 2*nodes; j++)
                {
                    sources[j]  = sgBA[j] + sgCO2[j] + 0.5*sc[j];
                }
                mZeroSource[celli] = sources[0];
                mOneSource[celli] = sources[1]; 
                mTwoSource[celli] = sources[2];
                mThreeSource[celli] = sources[3];
            }
            else
            {
                mZeroSource[celli] = ROOTVSMALL;
                mOneSource[celli] = ROOTVSMALL; 
                mTwoSource[celli] = ROOTVSMALL;
                mThreeSource[celli] = ROOTVSMALL;
            }   
        } //end of forAll
    }

    if 
    (
        PBEMethod == "QMOM" 
     && nNodes == 3 
    )
    {
        static int nodes = nNodes;
        int mOrder[2*nodes];
        for (int i = 0; i < 2*nNodes; i++){mOrder[i] = i;}
        double mom[2*nodes], we[nodes], vi[nodes], sgBA[2*nodes], 
        sgCO2[2*nodes], sc[2*nodes];
        double sources[2*nodes];

        double CO2_l_val, L_l_val, tmp_val, alphafoam;
        double LMAX, xBL_val, rhoPolyS;
        double cc1_val = 0.0;
        
        // Calling the model for density reaction mixture
        size_t T_pos = 
            modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
        size_t XOH_pos = 
            modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
        modena_model_argPos_check(density_reaction_mixturemodel);

        forAll(mesh.C(), celli)
        {
            CO2_l_val = wCO2_l[celli];
            L_l_val = wBA_l[celli];
            tmp_val = TS[celli];
            alphafoam = alpha2[celli];

            // Required variables for wCO2_Max
            double wCO2_g_val, wBA_g_val,prsr;
            wCO2_g_val= wCO2_g[celli]; 
            wBA_g_val = wBA_g[celli];
            prsr = p[celli];
            
            if (liquidMixtureDensitySurrogate)
            {
                // set input vector
                modena_inputs_set(inputs_den, Td_pos, TS[celli]);
                modena_inputs_set(inputs_den, XOHd_pos, XOH[celli]);
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
            
            if (alphafoam > 0.5)
            {
                mom[0] = mZero[celli];
                mom[1] = mOne[celli]; 
                mom[2] = mTwo[celli]; 
                mom[3] = mThree[celli]; 
                mom[4] = mFour[celli];  
                mom[5] = mFive[celli];

                // calling the functions, PDA, growthSource, coalescenceSource
                PDA(we, vi, mom, nodes);
                weight0[celli] = we[0];
                weight1[celli] = we[1];
                weight2[celli] = we[2];
                node0[celli] = vi[0];
                node1[celli] = vi[1]; 
                node2[celli] = vi[2]; 

                // Required variables for wCO2_Max
                double henry_Coefficient, bubble_radius, 
                   partialPressure_CO2, wCO2_Max;
                henry_Coefficient = henryCoefficient(tmp_val);
                bubble_radius = bubbleRadius(mom[0], mom[1]);
                partialPressure_CO2 = 
                    partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension, wCO2_g_val, wBA_g_val, 
                        prsr, bubble_radius
                    );
                wCO2_Max = 
                    wCO2Max
                    (
                        M_CO2, M_liq, partialPressure_CO2, henry_Coefficient
                    );

                if (blowingAgent == "n-pentane")
                {
                    LMAX = LliqMax (tmp_val);
                }
                if (blowingAgent == "R-11")
                {
                    xBL_val = xBL(tmp_val, dxdT);
                    LMAX = wBL_D(xBL_val, M_B, M_NCO, L0);
                }
                if (blowingAgent == "no")
                {
                    LMAX = 0.0;
                }
                growthSource
                (
                    sgBA, sgCO2, we, vi, nodes, mOrder, CO2_l_val, L_l_val, 
                    tmp_val, wCO2_Max, cc1_val, LMAX
                );
            
                g1_CO2[celli] = sgCO2[1];
                g1_BA[celli] = sgBA[1];
                cc1[celli] = cc1_val;
            
                if (XOH[celli] > XOH_Gel)
                {
                    for (int j = 0; j < 2*nodes; j++)
                    {
                        sc[j] = 0.0;
                    }
                }   
                else
                {
                    coalescenceSource(sc, we, vi, nodes, mOrder);
                }
            
                for (int j = 0; j < 2*nodes; j++)
                {
                    sources[j] = sgBA[j] + sgCO2[j] + 0.5*sc[j];
                }
                mZeroSource[celli] = sources[0];
                mOneSource[celli] = sources[1]; 
                mTwoSource[celli] = sources[2];
                mThreeSource[celli] = sources[3];
                mFourSource[celli] = sources[4];
                mFiveSource[celli] = sources[5];
            }
            else
            {
                mZeroSource[celli] = ROOTVSMALL;
                mOneSource[celli] = ROOTVSMALL; 
                mTwoSource[celli] = ROOTVSMALL;
                mThreeSource[celli] = ROOTVSMALL;
                mFourSource[celli] = ROOTVSMALL;
                mFiveSource[celli] = ROOTVSMALL;
            }   
        } //end of forAll
    }
}