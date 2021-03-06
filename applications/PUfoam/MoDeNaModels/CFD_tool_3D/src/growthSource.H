/**
@ingroup mod_3Dcfd
@file growthSource.H
@brief source term due to bubble growth
@fn void growthSource(double *sgBA, double *sgCO2, double *we, double *vi, int &nNodes, int *mOrder, double &CO2_l_val, double &L_l_val, double &tmp_val, double &wCO2_Max,double &cc1_val, double &LMAX, double *volumeGrowthBA, double *volumeGrowthCO2)
@param  sgBA source of growth due to blowing agent
@param  sgCO2 source of growth due to CO2
@param  we weights of quadrature
@param  vi nodes of quadrature
@param  nNodes number of nodes
@param  mOrder order of moments
@param  CO2_l_val weight fraction of CO2 in liquid
@param  L_l_val weight fraction of liquid blowing agent in liquid
@param  tmp_val temperature
@param  wCO2_Max weight fraction of maximum allowable CO2 in liquid
@param  cc1_val constant
@param  LMAX weight fraction of maximum allowable blowing agent in liquid
@return void
*/
void growthSource(double *, double *, double *, double *, int &, int *, double &, double &, double &, double &, double &, double &, double *, double *);

void growthSource
(
    double *sgBA, double *sgCO2, double *we, double *vi, int &nNodes, int *mOrder,
    double &CO2_l_val, double &L_l_val, double &tmp_val, double &wCO2_Max,
    double &cc1_val, double &LMAX,
    double *volumeGrowthBA, double *volumeGrowthCO2
)
{
    int i;
    int counter = 0;
    double k;

    if (L_l_val > LMAX)
    {
        cc1_val = 1.0;
    }
    else
    {
        cc1_val = 0.0;
    }

    while (counter < 2*nNodes)
    {
        sgBA[counter] = 0.0;
        sgCO2[counter] = 0.0;

        // using static_cast to convert moment order from int to double
        k = static_cast<double>(mOrder[counter]);

        if (counter == 0)
        {
            sgBA[counter] = 0.0;
            sgCO2[counter] = 0.0;
        }
        else if (counter == 1)
        {
            for (i = 0; i < nNodes; i++)
            {
                if(vi[i] > 0.0)
                {
                    sgBA[counter] += volumeGrowthBA[i]*we[i];
                    sgCO2[counter] += volumeGrowthCO2[i]*we[i];
                }
                else
                {
                    sgBA[counter] = sgBA[counter];
                    sgCO2[counter] = sgCO2[counter];
                }
            }
        }
        else
        {
            for (i = 0; i < nNodes; i++)
            {
                if (vi[i] > 0.0)
                {
                    sgBA[counter] += k*volumeGrowthBA[i]*we[i]*Foam::pow(vi[i], k - 1);
                    sgCO2[counter] += k*volumeGrowthCO2[i]*we[i]*Foam::pow(vi[i], k - 1);
                }
                else
                {
                    sgBA[counter] = sgBA[counter];
                    sgCO2[counter] = sgCO2[counter];
                }
            }
        }
        counter++;
    }
}