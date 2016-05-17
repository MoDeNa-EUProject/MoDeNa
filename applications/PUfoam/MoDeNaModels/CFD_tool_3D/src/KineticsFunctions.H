// Functions for kinetics calculations
//
// LliqMax - maximum allowable amount of liquid blowing agent (n-pentane) in liquid
// @param - double &tmptr, the input temperature.
// ddT_LliqMax - derivative of LliqMax with respect to temperature
// @param - double &tmptr, the input temperature.
// arrhenius - Arrhenius function
// @param - double& AOH - pre-exponential factor
// @param - double& EOH - activation energy
// @param - double& tempt - temperature
// BAdilution - dilution term for the blowing agent (n-pentane)
// @param - double& L_l - weight fraction of liquid blowing agent
// @param - double& rhoPoly - density of liquid mixture (polymer)
// @param - double& rhoBL - density of blowing agent
// QKinW - arrhenius term times by dilution term
// @param - similar to arrehenius and BAdilution functions
// henryCoefficient - Henry coefficient for CO2
// @param - double &T - Temperature
// bubbleRadius - radius of bubbles based on moments
// @param - const double m0 - moment of order zero
// @param - const double m1 - moment of order one
// partialPressureCO2 - partial pressure of CO2
// @param - double &M_CO2 - molecular weight of CO2
// @param - double &M_B - molecular weight of blowing agent
// @param - double &surfaceTension - surface tenstion
// @param - double &wCO2_g - weight fraction of CO2 in gas
// @param - double &wBA_g - weight frcation of blowing agent in gas
// @param - double &p - ambient pressure
// @param - double &R - bubble radius
// partialPressureBA - partial pressure of blowing agent
// @param - double &M_B - molecular weight of blowing agent
// @param - double &M_CO2 - molecular weight of CO2
// @param - double &surfaceTension - surface tenstion
// @param - double &wBA_g - weight frcation of blowing agent in gas
// @param - double &wCO2_g - weight fraction of CO2 in gas
// @param - double &p - ambient pressure
// @param - double &R - bubble radius

// wCO2Max - dissolved amount of CO2 in liquid
// @param - double &M_CO2 - molecular weight of CO2
// @param - double &M_liq - molecular weight of liquid mixture
// @param - double &pCO2 - partial pressure of CO2
// @param - double &henryCoeff - Henry coefficient
// xBL - mole fraction of blowing agent (R-11) in liquid polymer
// @param - double &T - temperature
// @param - double &dxdT - constant
// wBL_D - weight fraction of maximum allowable blowing agent (R-11) in liquid
// @param - double &xBL - mole fraction of blowing agent (R-11)
// @param - double &M_B - molecular weight of blowing agent (R-11)
// @param - double &M_NCO - molecular weight of NCO
// creamTemperature - The temperature that foaming process starts
// @param - double &xBL0 - initial mole fraction of the blowing agent (R-11)
// @param - double &dxdT - constant
// dTdtFirstOrder - first order accurate dTdt
// @param - Tc - current value of temperature
// @param - To - old value of temperature
// @param - timeStep - time step
// dTdtSecondOrder - second order accurate dTdt
// @param - Tc - current value of temperature
// @param - To - old value of temperature
// @param - Too - old-old value of temperature
// @param - timeStep - time step
// thermalDiffusivityGas - return thermal diffusivity of gas as a function of T
// @param = double &T - Temperature
// nodeRadius - bubble radius at the node
// @param - const double &v - bubble volume
double LliqMax (double &);
double ddT_LliqMax (double &);
double arrhenius (double &, double &, double &);
double BAdilution(double &, double &, double &);
double QKinW
    (
        double &, double &, double &, double &, double &, double &, double &
    );
double henryCoefficient(double &T);
double bubbleRadius (const double m0, const double m1);
double partialPressureCO2
    (
        double &M_CO2, double &M_B, double &surfaceTension, double &wCO2_g,
        double &wBA_g, double &p, double &R
    );
double partialPressureBA
    (
        double &M_B, double &M_CO2, double &surfaceTension,
        double &wBA_g_val, double &wCO2_g_val,
        double &prsr, double &bubble_radius
    );
double wCO2Max (double &M_CO2, double &M_liq, double &pCO2, double &henryCoeff);
double xBL (double &T, double &dxdT);
double wBL_D(double &xBL, double &M_B, double &M_NCO, double &L0);
double LmaxR11(double &tmptr);
double creamTemperature(double &xBL0, double &dxdT);
double dTdtFirstOrder(double &Tc, double &To, double &timeStep);
double dTdtSecondOrder(double &Tc, double &To, double &Too, double &timeStep);
double thermalDiffusivityGas(double &T);
double nodeRadius(const double &v);

double LliqMax (double& tmptr)
{
    // constants:
    const double a = 0.0064, h = 0.0551, T0 = 298.0, ww = 17.8;
    double lMax;
    if (tmptr > T0)
    {
        double tempDummy = Foam::pow((tmptr-T0),2.0);
        lMax = (a + h*Foam::exp((-tempDummy/(2.0*ww*ww))));
    }
    else
    {
        lMax = (a+h);
    }
    return (lMax);
}

double ddT_LliqMax (double& tmptr)
{
    // constants:
    const double h = 0.0551, T0 = 298.0, ww = 17.8;
    double tempDummy = Foam::pow((tmptr-T0),2.0);
    return ((-h*(tmptr-T0)*Foam::exp((-tempDummy/(2.0*ww*ww))))/(ww*ww));
}

double arrhenius (double& AOH, double& EOH, double& tempt)
{
    const double R = 8.3145; // J/mol K
    return (AOH*Foam::exp(-EOH/(R*tempt)));
}

double BAdilution(double& L_l, double& rhoPoly, double& rhoBL)
{
    return (1/(1+L_l*(rhoPoly/rhoBL)));
}

double QKinW
(
    double& AW, double& EW, double& tempi, double& Lliq,
    double& rhoPoly, double& rhoBL
)
{
    double R = 8.3145;
    return (AW*(Foam::exp(-EW/(R*tempi)))*(1/(1+Lliq*(rhoPoly/rhoBL))));
}

double henryCoefficient(double &T)
{
    if (T < 600)
    {
        // constants
        double a    = 1.771e7;
        double b    = -1.134e5;
        double c    = 320.2;
        double d    = -0.2563;

        return (a + b*T + c*T*T + d*T*T*T);
    }
    else
    {
        return (8.5e6);
    }
}

double bubbleRadius (const double m0, const double m1)
{
    if (m0 != 0.0 && m1 != 0.0)
    {
        double R;
        R   = Foam::pow((3.0*m1/(4.0*M_PI*m0)), 1.0/3.0);
        return R;
    }
    else
    {
        return (30e-6);
    }

}
double partialPressureCO2
(
    double &M_CO2, double &M_B, double &surfaceTension, double &wCO2_g,
    double &wBA_g, double &p, double &R
)
{
    double pCO2;
    if (wCO2_g == 0.0)
    {
        pCO2 = SMALL;
    }
    else
    {
        pCO2 =
        (
            (wCO2_g/M_CO2)/(wBA_g/M_B + wCO2_g/M_CO2)
        );
    }

    return (1e5*pCO2);
}
double partialPressureBA
(
    double &M_B, double &M_CO2, double &surfaceTension, double &wBA_g,
    double &wCO2_g, double &p, double &R
)
{
    double pBA;
    if (wBA_g == 0.0)
    {
        pBA = SMALL;
    }
    else
    {
        pBA =

            (
                (wBA_g/M_B)/(wBA_g/M_B + wCO2_g/M_CO2)
            );
    }
    return (1e5*pBA);
}

double wCO2Max (double &M_CO2, double &M_liq, double &pCO2, double &henryCoeff)
{
    if ((henryCoeff - pCO2 ) > 0.0)
    {
        return ((M_CO2/M_liq)*(pCO2/(henryCoeff - pCO2)));
    }
    else
    {
        Info<< "\nWarning! Invalid wCO2Max value!" << endl;
        Info<< "'wCO2Max' is replaced by a constant." << endl;
        return (4.4e-4);
    }

}

double xBL (double &T, double &dxdT)
{
    double xBL_value;
    xBL_value   = dxdT*(T - 300.0) + 0.5;

    if (xBL_value < 0.0)
    {
        return (0.0);
    }
    else if (xBL_value > 0.5)
    {
        return (0.5);
    }
    else
    {
        return (xBL_value);
    }
}

double wBL_D(double &xBL, double &M_B, double &M_NCO, double &L0)
{
    double Lm;
    Lm  = (xBL/(1.0 - xBL))*(M_B/M_NCO);

    if (Lm < 1.0e-4)
    {
        Lm = 1.0e-4;
    }
    return (Lm);
}

double LmaxR11(double &tmptr)
{
    // constants:
        const double a = 1e-7, h = 4.293, T0 = 203.4, ww = 40.02;
        double lMax;
        if (tmptr > T0)
        {
            double tempDummy = Foam::pow((tmptr-T0),2.0);
            lMax = (a + h*Foam::exp((-tempDummy/(2.0*ww*ww))));
        }
        else
        {
            lMax = (a+h);
        }
    return (lMax);
}

double creamTemperature(double &xBL0, double &dxdT)
{
    return ((xBL0-0.5)/dxdT + 300.0);
}

double dTdtFirstOrder(double &Tc, double &To, double &timeStep)
{
    return (mag(Tc - To)/(max(timeStep,ROOTVSMALL)));
}
double dTdtSecondOrder(double &Tc, double &To, double &Too, double &timeStep)
{
    return (mag((3.0*Tc - 4*To + Too))/(max((2.0*timeStep),ROOTVSMALL)));
}

double thermalDiffusivityGas(double &T)
{
    return (2.56008e-3*std::pow(T,0.5)/(1+207.009/T))/1010/1.2;
}
double nodeRadius(const double &v)
{
    return (Foam::pow(3.0*v/(4.0*M_PI), 1/3.0));
}
