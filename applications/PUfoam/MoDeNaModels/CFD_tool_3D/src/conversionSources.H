scalar temp_current, 
       temp_old, 
       Arhns_current, 
       Arhns_old, 
       rhoFoam_current, 
       rhoFoam_old, 
       alphaFoam_current, 
       alphaFoam_old, 
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
    rhoFoam_current = rho_foam[celli];
    rhoFoam_old = rho_foam.oldTime()[celli];
    alphaFoam_current = alpha2[celli];
    alphaFoam_old = alpha2.oldTime()[celli];

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