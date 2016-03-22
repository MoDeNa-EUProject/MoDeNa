if (KineticsModel == "generic")
{
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
}
if (KineticsModel == "simple")
{
    // Calling the simpleKinetics model and retrieving the source terms
    size_t EG_NCO_Pos = 
        modena_model_inputs_argPos(kinetics, "'EG_NCO'");
    size_t EG_OH_Pos = 
        modena_model_inputs_argPos(kinetics, "'EG_OH'");
    size_t H2O_Pos = 
        modena_model_inputs_argPos(kinetics, "'H2O'");
    size_t CO2_Pos = 
        modena_model_inputs_argPos(kinetics, "'CO2'");
    size_t PENTANE_Pos = 
        modena_model_inputs_argPos(kinetics, "'PENTANE'");
    size_t POLYMER_Pos = 
        modena_model_inputs_argPos(kinetics, "'POLYMER'");
    size_t POLYMERBLOW_Pos = 
        modena_model_inputs_argPos(kinetics, "'POLMERBLOW'");
    size_t UREA_Pos = 
        modena_model_inputs_argPos(kinetics, "'UREA'");
    size_t R_1_temp_Pos = 
        modena_model_inputs_argPos(kinetics, "'R_1_temp'");
    // outputs argPos
    size_t source_EG_NCO_Pos = 
        modena_model_outputs_argPos(kinetics, "source_EG_NCO");
    size_t source_EG_OH_Pos = 
        modena_model_outputs_argPos(kinetics, "source_EG_OH");
    size_t source_H2O_Pos = 
        modena_model_outputs_argPos(kinetics, "source_H2O");
    size_t source_CO2_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CO2");
    size_t source_PENTANE_Pos = 
        modena_model_outputs_argPos(kinetics, "source_PENTANE");
    size_t source_POLYMER_Pos = 
        modena_model_outputs_argPos(kinetics, "source_POLYMER");
    size_t source_POLYMERBLOW_Pos = 
        modena_model_outputs_argPos(kinetics, "source_POLMERBLOW");
    size_t source_UREA_Pos = 
        modena_model_outputs_argPos(kinetics, "source_UREA");
    size_t source_R_1_temp_Pos = 
        modena_model_outputs_argPos(kinetics, "source_R_1_temp");
    modena_model_argPos_check(kinetics);

    forAll(mesh.C(), celli)
    {
        double alphaFoam, rhoFoam;
        alphaFoam = alpha2[celli];
        rhoFoam = rho_foam[celli];
        
        // set input vector
        modena_inputs_set(inputs_kinetics, EG_NCO_Pos, EG_NCO[celli]);
        modena_inputs_set(inputs_kinetics, EG_OH_Pos, EG_OH[celli]);
        modena_inputs_set(inputs_kinetics, H2O_Pos, H2O[celli]);
        modena_inputs_set(inputs_kinetics, CO2_Pos, CO2[celli]);
        modena_inputs_set(inputs_kinetics, PENTANE_Pos, PENTANE[celli]);
        modena_inputs_set(inputs_kinetics, POLYMER_Pos, POLYMER[celli]);
        modena_inputs_set(inputs_kinetics, POLYMERBLOW_Pos, POLMERBLOW[celli]);
        modena_inputs_set(inputs_kinetics, UREA_Pos, UREA[celli]);
        modena_inputs_set(inputs_kinetics, R_1_temp_Pos, R_1_temp[celli]);
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
        // source terms for the simpleKinetics surrogate model
        source_EG_NCO[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_EG_NCO_Pos);
        source_EG_NCO[celli] = 
            (source_EG_NCO[celli] < 0.0) ? 0.0 : source_EG_NCO[celli];
        source_EG_OH[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_EG_OH_Pos);
        source_EG_OH[celli] = 
            (source_EG_OH[celli] < 0.0) ? 0.0 : source_EG_OH[celli];
        source_H2O[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_H2O_Pos);
        source_H2O[celli] = 
            (source_H2O[celli] < 0.0) ? 0.0 : source_H2O[celli];
        source_CO2[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_CO2_Pos);
        source_CO2[celli] = 
            (source_CO2[celli] < 0.0) ? 0.0 : source_CO2[celli];
        source_PENTANE[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_PENTANE_Pos);
        source_PENTANE[celli] = 
            (source_PENTANE[celli] < 0.0) ? 0.0 : source_PENTANE[celli];
        source_POLYMER[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_POLYMER_Pos);
        source_POLYMER[celli] = 
            (source_POLYMER[celli] < 0.0) ? 0.0 : source_POLYMER[celli];
        source_POLMERBLOW[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_POLYMERBLOW_Pos);
        source_POLMERBLOW[celli] = 
            (source_POLMERBLOW[celli] < 0.0) ? 0.0 : source_POLMERBLOW[celli];
        source_UREA[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_UREA_Pos);
        source_UREA[celli] = 
            (source_UREA[celli] < 0.0) ? 0.0 : source_UREA[celli];
        source_R_1_temp[celli] = 
            rhoFoam*alphaFoam*modena_outputs_get(outputs_kinetics, source_R_1_temp_Pos);
        source_R_1_temp[celli] = 
            (source_R_1_temp[celli] < 0.0) ? 0.0 : source_R_1_temp[celli];
    }
}
if (KineticsModel == "RF-1")
{
    // inputs argPos
    size_t Catalyst_1_Pos = modena_model_inputs_argPos(kinetics, "'Catalyst_1'");
    size_t CE_A0_Pos = modena_model_inputs_argPos(kinetics, "'CE_A0'");
    size_t CE_A1_Pos = modena_model_inputs_argPos(kinetics, "'CE_A1'");
    size_t CE_B_Pos = modena_model_inputs_argPos(kinetics, "'CE_B'");
    size_t CE_B2_Pos = modena_model_inputs_argPos(kinetics, "'CE_B2'");
    size_t CE_I0_Pos = modena_model_inputs_argPos(kinetics, "'CE_I0'");
    size_t CE_I1_Pos = modena_model_inputs_argPos(kinetics, "'CE_I1'");
    size_t CE_I2_Pos = modena_model_inputs_argPos(kinetics, "'CE_I2'");
    size_t CE_PBA_Pos = modena_model_inputs_argPos(kinetics, "'CE_PBA'");
    size_t CE_Breac_Pos = modena_model_inputs_argPos(kinetics, "'CE_Breac'");
    size_t CE_Areac0_Pos = modena_model_inputs_argPos(kinetics, "'CE_Areac0'");
    size_t CE_Areac1_Pos = modena_model_inputs_argPos(kinetics, "'CE_Areac1'");
    size_t CE_Ireac0_Pos = modena_model_inputs_argPos(kinetics, "'CE_Ireac0'");
    size_t CE_Ireac1_Pos = modena_model_inputs_argPos(kinetics, "'CE_Ireac1'");
    size_t CE_Ireac2_Pos = modena_model_inputs_argPos(kinetics, "'CE_Ireac2'");
    size_t Bulk_Pos = modena_model_inputs_argPos(kinetics, "'Bulk'");
    size_t R_1_Pos = modena_model_inputs_argPos(kinetics, "'R_1'");
    size_t R_1_mass_Pos = modena_model_inputs_argPos(kinetics, "'R_1_mass'");
    size_t R_1_temp_RF1_Pos = modena_model_inputs_argPos(kinetics, "'R_1_temp'");
    size_t R_1_vol_RF1_Pos = modena_model_inputs_argPos(kinetics, "'R_1_vol'");

    // outputs argPos
    size_t source_Catalyst_1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_Catalyst_1");
    size_t source_CE_A0_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_A0");
    size_t source_CE_A1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_A1");
    size_t source_CE_B_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_B");
    size_t source_CE_B2_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_B2");
    size_t source_CE_I0_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_I0");
    size_t source_CE_I1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_I1");
    size_t source_CE_I2_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_I2");
    size_t source_CE_PBA_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_PBA");        
    size_t source_CE_Breac_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_Breac");
    size_t source_CE_Areac0_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_Areac0");
    size_t source_CE_Areac1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_Areac1");
    size_t source_CE_Ireac0_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_Ireac0");
    size_t source_CE_Ireac1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_Ireac1");
    size_t source_CE_Ireac2_Pos = 
        modena_model_outputs_argPos(kinetics, "source_CE_Ireac2");
    size_t source_Bulk_Pos = 
        modena_model_outputs_argPos(kinetics, "source_Bulk");
    size_t source_R_1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_R_1");
    size_t source_R_1_mass_Pos = 
        modena_model_outputs_argPos(kinetics, "source_R_1_mass");
    size_t source_R_1_temp_RF1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_R_1_temp");
    size_t source_R_1_vol_RF1_Pos = 
        modena_model_outputs_argPos(kinetics, "source_R_1_vol");

    modena_model_argPos_check(kinetics);
    
    forAll(mesh.C(), celli)
    {
        double alphaFoam, rhoFoam;
        alphaFoam = alpha2[celli];
        rhoFoam = rho_foam[celli];
        // set input vector

        modena_inputs_set(inputs_kinetics, Catalyst_1_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*Catalyst_1[celli]));
        modena_inputs_set(inputs_kinetics, CE_A0_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_A0[celli]));
        modena_inputs_set(inputs_kinetics, CE_A1_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_A1[celli]));
        modena_inputs_set(inputs_kinetics, CE_B_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_B[celli]));
        modena_inputs_set(inputs_kinetics, CE_B2_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_B2[celli]));
        modena_inputs_set(inputs_kinetics, CE_I0_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_I0[celli]));
        modena_inputs_set(inputs_kinetics, CE_I1_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_I1[celli]));
        modena_inputs_set(inputs_kinetics, CE_I2_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_I2[celli]));
        modena_inputs_set(inputs_kinetics, CE_PBA_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_PBA[celli]));
        modena_inputs_set(inputs_kinetics, CE_Breac_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_Breac[celli]));
        modena_inputs_set(inputs_kinetics, CE_Areac0_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_Areac0[celli]));
        modena_inputs_set(inputs_kinetics, CE_Areac1_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_Areac1[celli]));
        modena_inputs_set(inputs_kinetics, CE_Ireac0_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_Ireac0[celli]));
        modena_inputs_set(inputs_kinetics, CE_Ireac1_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_Ireac1[celli]));
        modena_inputs_set(inputs_kinetics, CE_Ireac2_Pos, (rhoPoly*(1/scalar(1000))*(1/2.65000e-1)*CE_Ireac2[celli]));
        modena_inputs_set(inputs_kinetics, Bulk_Pos, Bulk[celli]);
        modena_inputs_set(inputs_kinetics, R_1_Pos, R_1[celli]);
        modena_inputs_set(inputs_kinetics, R_1_mass_Pos, R_1_mass[celli]);
        modena_inputs_set(inputs_kinetics, R_1_temp_RF1_Pos, (R_1_temp[celli]-273.15));
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
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_Catalyst_1_Pos);
        // source_Catalyst_1[celli] = 
        //     (source_Catalyst_1[celli] < 0.0) ? 0.0 : source_Catalyst_1[celli];
        Info<< "alphaFoam = " << alpha2[celli] << endl;
        source_CE_A0[celli] = 
            1000.0*alphaFoam*2.65000e-1
            *modena_outputs_get(outputs_kinetics, source_CE_A0_Pos);
        Info<< "source_CE_A0[celli] with multiplication: " << source_CE_A0[celli] << endl;
        Info<< "source_CE_A0[celli] without multiplication: " << modena_outputs_get(outputs_kinetics, source_CE_A0_Pos) << endl; 
        // source_CE_A0[celli] =
        //     (source_CE_A0[celli] < 0.0) ? 0.0 : source_CE_A0[celli];
        
        source_CE_A1[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_A1_Pos);
        // source_CE_A1[celli] =
        //     (source_CE_A1[celli] < 0.0) ? 0.0 : source_CE_A1[celli];

        source_CE_B[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_B_Pos);
        // source_CE_B[celli] =
        //     (source_CE_B[celli] < 0.0) ? 0.0 : source_CE_B[celli];

        source_CE_B2[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_B2_Pos);
        // source_CE_B2[celli] =
        //     (source_CE_B2[celli] < 0.0) ? 0.0 : source_CE_B2[celli];

        source_CE_I0[celli] = 
        scalar(1000)*alphaFoam*scalar(2.65000e-1)
        *modena_outputs_get(outputs_kinetics, source_CE_I0_Pos);
        // source_CE_I0[celli] =
        //     (source_CE_I0[celli] < 0.0) ? 0.0 : source_CE_I0[celli];
        
        source_CE_I1[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_I1_Pos);
        // source_CE_I1[celli] =
        //     (source_CE_I1[celli] < 0.0) ? 0.0 : source_CE_I1[celli];

        source_CE_I2[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_I2_Pos);
        // source_CE_I2[celli] =
        //     (source_CE_I2[celli] < 0.0) ? 0.0 : source_CE_I2[celli];

        source_CE_PBA[celli] = 
            scalar(1000)*scalar(2.65000e-1)*alphaFoam
            *modena_outputs_get(outputs_kinetics, source_CE_PBA_Pos);
        // source_CE_PBA[celli] =
        //     (source_CE_PBA[celli] < 0.0) ? 0.0 : source_CE_PBA[celli];

        source_CE_Breac[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_Breac_Pos);
        // source_CE_Breac[celli] =
        //     (source_CE_Breac[celli] < 0.0) ? 0.0 : source_CE_Breac[celli];

        source_CE_Areac0[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_Areac0_Pos);
        // source_CE_Areac0[celli] =
        //     (source_CE_Areac0[celli] < 0.0) ? 0.0 : source_CE_Areac0[celli];

        source_CE_Areac1[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_Areac1_Pos);
        // source_CE_Areac1[celli] =
        //     (source_CE_Areac1[celli] < 0.0) ? 0.0 : source_CE_Areac1[celli];

        source_CE_Ireac0[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_Ireac0_Pos);
        // source_CE_Ireac0[celli] =
        //     (source_CE_Ireac0[celli] < 0.0) ? 0.0 : source_CE_Ireac0[celli];

        source_CE_Ireac1[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_Ireac1_Pos);
        // source_CE_Ireac1[celli] =
        //     (source_CE_Ireac1[celli] < 0.0) ? 0.0 : source_CE_Ireac1[celli];

        source_CE_Ireac2[celli] = 
            scalar(1000)*alphaFoam*scalar(2.65000e-1)
            *modena_outputs_get(outputs_kinetics, source_CE_Ireac2_Pos);
        // source_CE_Ireac2[celli] =
        //     (source_CE_Ireac2[celli] < 0.0) ? 0.0 : source_CE_Ireac2[celli];

        source_Bulk[celli] = 
            modena_outputs_get(outputs_kinetics, source_Bulk_Pos);
        // source_Bulk[celli] =
        //     (source_Bulk[celli] < 0.0) ? 0.0 : source_Bulk[celli];

        source_R_1[celli] = 
            modena_outputs_get(outputs_kinetics, source_R_1_Pos);
        // source_R_1[celli] =
        //     (source_R_1[celli] < 0.0) ? 0.0 : source_R_1[celli];

        source_R_1_mass[celli] = 
            modena_outputs_get(outputs_kinetics, source_R_1_mass_Pos);
        // source_R_1_mass[celli] =
        //     (source_R_1_mass[celli] < 0.0) ? 0.0 : source_R_1_mass[celli];

        source_R_1_temp_RF1[celli] = 
            modena_outputs_get(outputs_kinetics, source_R_1_temp_RF1_Pos);
        // source_R_1_temp_RF1[celli] =
        //     (source_R_1_temp_RF1[celli] < 0.0) ? 0.0 : source_R_1_temp_RF1[celli];

        source_R_1_vol_RF1[celli] = 
            modena_outputs_get(outputs_kinetics, source_R_1_vol_RF1_Pos);
        // source_R_1_vol_RF1[celli] =
        //     (source_R_1_vol_RF1[celli] < 0.0) ? 0.0 : source_R_1_vol_RF1[celli];
    
    // Info<< "source_Catalyst_1[celli] = "    << source_Catalyst_1[celli] << endl;
    // Info<< "source_CE_A0[celli] = "         << source_CE_A0[celli] << endl;
    // Info<< "source_CE_A1[celli] = "         << source_CE_A1[celli] << endl;
    // Info<< "source_CE_B[celli] = "          << source_CE_B[celli] << endl;
    // Info<< "source_CE_B2[celli] = "         << source_Catalyst_1[celli] << endl;
    // Info<< "ssource_CE_I0[celli] = "        << source_CE_I0[celli] << endl;
    // Info<< "source_CE_I0[celli] = "         << source_CE_I0[celli] << endl;
    // Info<< "source_CE_I1[celli] = "         << source_CE_I1[celli] << endl;
    // Info<< "source_CE_I2[celli] = "         << source_CE_I2[celli] << endl;
    // Info<< "source_CE_PBA[celli] = "        << source_CE_PBA[celli] << endl;
    // Info<< "source_CE_Breac[celli] = "      << source_CE_Breac[celli] << endl;
    // Info<< "source_CE_Areac0[celli] = "     << source_CE_Areac0[celli] << endl;
    // Info<< "source_CE_Areac1[celli] = "     << source_CE_Areac1[celli] << endl;
    // Info<< "source_CE_Ireac0[celli] = "     << source_CE_Ireac0[celli] << endl;
    // Info<< "source_CE_Ireac1[celli] = "     << source_CE_Ireac1[celli] << endl;
    // Info<< "source_CE_Ireac2[celli] = "     << source_CE_Ireac2[celli] << endl;
    // Info<< "source_Bulk[celli] = "          << source_Bulk[celli] << endl;
    // Info<< "source_R_1[celli] = "           << source_R_1[celli] << endl;
    // Info<< "source_R_1_mass[celli] = "      << source_R_1_mass[celli] << endl;
    // Info<< "source_R_1_vol_RF1[celli] = "   << source_R_1_vol_RF1[celli] << endl;
    Info<< endl;
    }

}