double dXOHdt,XOH_C,XOH_Old,dWdt,XW_C,XW_Old,dLdt,alphaFO,rhoFO,
    C_TOT, L_l_C, L_l_Old, T_C,T_O, wBA_g_val, xBL_val; 

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
    
    T_C = mixture.T()[celli]; 
    T_O = mixture.T().oldTime()[celli]; 

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
        lmax = LliqMax(T_C);
        if (lmax < L0)
        {
            dLdt = 
                -1.0*ddT_LliqMax(T_C)
                *((T_C-T_O)/(max(runTime.deltaTValue(),ROOTVSMALL)));
        }
        else
        {
            dLdt = 0.0;
        }   
    }

    if (blowingAgent == "n-pentane")
    {
        C_TOT = C_Poly;
        if 
        (
            alphaFO > 0.5
         && XOH[celli] < XOH_Gel
        ) 
        {
            TSource[celli] = 
            (
                (alphaFO*rhoFO)*((-DH_OH*COH_0*dXOHdt)/(rhoPoly*C_TOT)
                + (-DH_W*CW_0*dWdt)/(rhoPoly*C_TOT) + (latenth*dLdt)/(rhoPoly*C_TOT))
            );
        }
        else if 
        (
            alphaFO > 0.5
         && XOH[celli] >= XOH_Gel
        )
        {
            TSource[celli] = (alphaFO*rhoFO)*(latenth*dLdt)/(rhoPoly*C_TOT);
        }
        else
        {
            TSource[celli] = ROOTVSMALL;
        }   
    }   
    
    if (blowingAgent == "R-11")
    {
        if (wBA_g_val > 0.0)
        {
            xBL_val = xBL(T_C, dxdT);
            
            C_TOT = 
            (
                C_Poly + (-(M_B/M_NCO)*(1.0/
                (Foam::pow((1.0 - xBL_val),2)))*(dxdT))*latenth
            );
            
            if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
            {
                TSource[celli] = 
                (
                    (alphaFO*rhoFO)*((-DH_OH*COH_0*dXOHdt)/(rhoPoly*C_TOT)
                  + (-DH_W*CW_0*dWdt)/(rhoPoly*C_TOT))
                );
            }
            else if (alphaFO > 0.5 && XOH[celli] >= XOH_Gel)
            {
                TSource[celli] = 
                (
                    (alphaFO*rhoFO)*((-(M_B/M_NCO)*
                    (1.0/(Foam::pow((1.0 - xBL_val),2)))*(dxdT))*latenth)
                );
            }
            else
            {
                TSource[celli] = ROOTVSMALL;
            }
        }
        else
        {
            C_TOT = C_Poly;
            if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
            {
                TSource[celli] = 
                (
                    (alphaFO*rhoFO)*((-DH_OH*COH_0*dXOHdt)/
                    (rhoPoly*C_TOT) + (-DH_W*CW_0*dWdt)/(rhoPoly*C_TOT))
                );
            }
            else
            {
                TSource[celli] = ROOTVSMALL;
            }
        }
    }
    if (blowingAgent == "no")
    {
        C_TOT = C_Poly;
        if (alphaFO > 0.5 && XOH[celli] < XOH_Gel)
        {
            TSource[celli] = 
            (
                (alphaFO*rhoFO)*((-DH_OH*COH_0*dXOHdt)/
                (rhoPoly*C_TOT) + (-DH_W*CW_0*dWdt)/(rhoPoly*C_TOT) )
            );
        }
        else
        {
            TSource[celli] = ROOTVSMALL;
        }   
    }   
    if (blowingAgent == "n-pentane")
    {
        if (alphaFO > 0.5)
        {   
            if (rhoFO > 48)
            {
                thermalConductivity[celli] = 
                (
                    (scalar(8.7006e-8)*Foam::pow(rhoFO,2) + 
                    scalar(8.4674e-5)*rhoFO+scalar(1.16e-2))
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
            thermalConductivity[celli] = 0.03;
            thermalDiff_foam[celli] = (thermalConductivity[celli]/(rhoFO*C_Poly));
        }
        else
        {
            thermalDiff_foam[celli] = thermalDiffusivityGas(T_C);
        }
        thermalDiff_gas[celli] = thermalDiffusivityGas(T_C);
    }
        thermalDiffusivity[celli] = 
        (
            alpha1[celli]*thermalDiff_gas[celli] 
          + alpha2[celli]*thermalDiff_foam[celli];
        );
}
