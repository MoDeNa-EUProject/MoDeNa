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
       QW;

forAll(mesh.C(), celli)
{
    temp_current = TS[celli];
    temp_old = TS.oldTime()[celli];  
    rhoFoam_current = rho_foam[celli];
    rhoFoam_old = rho_foam.oldTime()[celli];
    alphaFoam_current = alpha2[celli];
    alphaFoam_old = alpha2.oldTime()[celli];

    if (XW[celli] < 0.0)
    {
        XW[celli] = ROOTVSMALL;
    }
    if (XW[celli] > 1.0)
    {
        XW[celli] = ROOTVSMALL;
    }
    QW = QKinW (AW, EW, temp_current, L0, rhoPoly, rhoBL);
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
        BA_dilution = BAdilution(L0,rhoPoly,rhoBL);
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