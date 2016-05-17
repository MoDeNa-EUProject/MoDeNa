/** @file partialPressure.h
    @brief functions related to the calculations of partial pressures
    @fn void ddtpartialPressure(const state_type &y , const double t , const double dt , double *dpdt , double *pOld, const double p_1, const double p_2, const double R)
    @brief time derivative of partial pressure
    @param  const state_type &y -  vector of all the variables
    @param  const double t - time
    @param  double *dpdt - time derivative of pressure
    @param  double *pOld - pressure values at the previous time step
    @param  const double p_1 - partial pressure of blowing agent   
    @param  const double p_2 - partial pressure of CO2
    @param  const double R - bubble radius
    @return void
    @fn double partialPressureBA(const state_type &y)
    @brief partial pressure of the physical blowing agent
    @param const state_type &y vector of all the variables
    @return partial pressure of the physical blowing agent
    @fn double partialPressureCO2(const state_type &y)
    @brief partial pressure of CO2
    @param const state_type &y - vector of all the variables 
    @return partial pressure of CO2
*/
void 	ddtpartialPressure(const state_type &y , const double t , const double dt , double *dpdt , double *pOld, const double p_1, const double p_2, const double R);
bool isNaN(double var);
double bubbleRadius (const double m0, const double m1);
double nodeRadius(const double &v);
double 	partialPressureBA(const state_type &y);
double	partialPressureCO2(const state_type &y);

void ddtpartialPressure(const state_type &y , const double t , const double dt , double *dpdt , double *pOld, const double p_1, const double p_2, const double R)
{
	// dp1dt:
    dpdt[0] = (p_1 - pOld[0])/(dt);
	pOld[0] = p_1;

    // dp2dt:
    dpdt[1] = (p_2 - pOld[1])/(dt);
	pOld[1] = p_2;
}
bool isNaN(double var)
{
    volatile double d = var;
    return d != d;
}
double bubbleRadius (const double m0, const double m1)
{
    double R;
    R   = pow((3.0*m1/(4.0*M_PI*m0)), 1.0/3.0);
    if (isNaN(R)) {
        R=init_size;
    }
    return R;
}

double nodeRadius(const double &v)
{
    return (pow((6.0*v/M_PI),1.0/3)/2.0);
}
double partialPressureBA(const state_type &y)
{
    double L_g      = y[4];
    double CO2_g    = y[6];
    double m0       = y[7];
    double m1       = y[8];
    // double PENTANE  = y[15];
    // double CO2      = y[14];

    double R = bubbleRadius(m0, m1);

    double p_1;
    if (L_g == 0.0)
    {
        p_1     = 0.0;
    }
    else
    {
        p_1     = ((L_g/M_B)/(L_g/M_B + CO2_g/M_CO2 + air_g/M_air)) * (Pr + 2*surfaceTension/R);
    }

    return p_1;
}

double partialPressureCO2(const state_type &y)
{
    double L_g      = y[4];
    double CO2_g    = y[6];
    double m0       = y[7];
    double m1       = y[8];
    // double PENTANE  = y[15];
    // double CO2      = y[14];

    double R = bubbleRadius(m0, m1);

    double p_2;

    if (CO2_g == 0.0)
    {
        p_2     = 0.0;
    }
    else
    {
        p_2     = ((CO2_g/M_CO2)/(L_g/M_B + CO2_g/M_CO2 + air_g/M_air)) * (Pr + 2*surfaceTension/R);
    }

    return p_2;
}