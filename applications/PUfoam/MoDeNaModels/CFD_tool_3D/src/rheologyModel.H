// Constant viscosity (Newtonian)
if (viscosityModel == "constant")
{
    forAll(mesh.C(), celli)
    {
        muFoamCorr[celli] = 1.0e-3;
    }
}

// Castro-Macosko model (Newtonian)
if (viscosityModel == "castro-macosko")
{
    double muTemp,alphacell,CXOH,CXW;
    double muinf = 10.3*1e-8;
    double XNCOGel = 1.0;
    double EnuR = 4970.0;
    double nua = 1.5;
    double nub = 1.0;
    double nuc = 0.0;

    forAll(mesh.C(), celli)
    {
        muTemp = TS[celli];
        alphacell = alpha2[celli];
        CXOH = XOH[celli];
        CXW = XW[celli];

        if (alphacell > 0.5)
        {
            XNCO[celli] = (CNCO_0 - (CNCO_0 - CXOH*COH_0 - 2*CXW*CW_0))/CNCO_0;
            muFoamCorr[celli] =
                (
                (muinf*Foam::exp(EnuR/(muTemp)) *
                Foam::pow((XNCOGel/max((XNCOGel - XNCO[celli]),ROOTVSMALL)),
                (nua + nub*XNCO[celli] + nuc*XNCO[celli]*XNCO[celli])))
                );
        }
    }
}

// Modified Bird-Carreau model (non-Newtonian)
if (viscosityModel == "bird-carreau")
{

    double coeffalpha,coeffn,muTemp,lambdaBC,XOHC,alphacell,a,b,c,d,mu0l,muinfl;
    scalar Amu = 0.0387;
    scalar EmuR = 10000.0;
    // model constant
    a = 1.5;
    b = 1.0;
    c = 0.0;
    d = 0.001;
    mu0l = 0.195;
    muinfl = 0.266;
    coeffalpha = 2.0;
    lambdaBC = 11.35;
    coeffn = 0.2;

    forAll(mesh.C(), celli)
    {
        alphacell = alpha2[celli];
        if (alphacell > 0.5)
        {
            muTemp = TS[celli];
            XOHC = XOH[celli];

            volScalarField shearRate = Foam::mag(fvc::grad(U));

            mu0[celli] =
            (
                (Foam::log(XOHC+d) - Foam::log(d) +
                Foam::pow((XOH_Gel/max((XOH_Gel - XOHC),0.00001)),
                ((a+b*XOHC+c*XOHC*XOHC))))*mu0l
            );
            muinf[celli] =
            (
                (Foam::log(XOHC + d) - Foam::log(d) +
                Foam::pow((XOH_Gel/max((XOH_Gel - XOHC),0.00001)),
                ((a + b*XOHC + c*XOHC*XOHC))))*muinfl
            );

            muFoam[celli] =
            (
                ((muinf[celli] + (mag((mu0[celli] - muinf[celli]))*
                Foam::pow((1.0 + Foam::pow((shearRate[celli]*lambdaBC),coeffalpha)),
                ((coeffn - 1)/coeffalpha)))))
            );
            muFoamCorr[celli] = (muFoam[celli]*(Amu*Foam::exp(EmuR/RR/muTemp)));
        }
        else
        {
            muFoamCorr[celli] = 1.0e-5;
        }
    }
}

if (viscosityModel == "MoDeNaRheology")
{
    forAll(mesh.C(), celli)
    {
        if (alpha2[celli] > 0.5)
        {
            volScalarField shearRate = Foam::mag(fvc::grad(U));
                // set input vector
            modena_inputs_set(inputs_rheo, temp_rheopos, TS[celli]);
            modena_inputs_set(inputs_rheo, conv_rheopos, XOH[celli]);
            modena_inputs_set(inputs_rheo, shear_rheopos, shearRate[celli]);
            modena_inputs_set(inputs_rheo, m0_rheopos, mZero[celli]);
            modena_inputs_set(inputs_rheo, m1_rheopos, mOne[celli]);
            // call the model
            int ret_rheo = modena_model_call (rheologymodel, inputs_rheo, outputs_rheo);
            if (modena_error_occurred())
            {
                modena_inputs_destroy (inputs_rheo);
                modena_outputs_destroy (outputs_rheo);
                modena_model_destroy (rheologymodel);
                Info<< "MoDeNa Error: " << (modena_error()) << endl;
                // exit(modena_error());
            }
            muFoamCorr[celli] = modena_outputs_get(outputs_rheo, 0);
        }
    }
}
muMixture = alpha1*muAir + alpha2*muFoamCorr;
muMixture.correctBoundaryConditions();

