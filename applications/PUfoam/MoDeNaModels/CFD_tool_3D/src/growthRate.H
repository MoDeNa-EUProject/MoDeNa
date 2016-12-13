// three different functions for growth rate
// - growthRateConst() : constant growth rate
// - growthRateCO2() : growth rate due to CO2
// @param - double &CO2_l - weight fraction of CO2 in liquid
// @param - double &wCO2_Max - weight fraction of maximum allowable CO2 in liquid
// - growthRateBA() : growth rate due to blowing agent
// @param - double &L_l - weight fraction of liquid blowing agent in liquid
// @param - double &LMAX - weight fraction of maximum allowable blowing agent in liquid
double growthRateConst();
double growthRateCO2(double &, double &);
double growthRateBA (double &, double &);

double growthRateConst()
{
    return 0.0955;
}

double growthRateCO2(double &CO2_l, double &wCO2_Max)
{
    double G0 = 1.0e-14;
    if (CO2_l > wCO2_Max)
    {
        return (G0*(CO2_l - wCO2_Max)/max(wCO2_Max,4.4e-4));
    }
    else
    {
        return (0.0);
    }
}

double growthRateBA(double &L_l, double &LMAX)
{
    double G0 = 1.0e-14;

    if (LMAX != 0.0 && L_l > LMAX)
    {
        return (G0*(L_l - LMAX)/LMAX);
    }
    else
    {
        return (0.0);
    }
}
