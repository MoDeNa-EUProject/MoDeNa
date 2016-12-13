/** @file growth.h
	@brief Source terms due to the bubble growth
	@fn void growthSource(double *sgBA, double *sgCO2, double *we, double *vi, int &nNodes, int *mOrder, double &CO2_l, double &L_l, double &tmp, double &dVdt_2, double &dVdt_1)
	@param  double *sgBA - source terms due to blowing agent
	@param  double *sgCO2 - source terms due to CO2
	@param  double *we - weights of quadrature approximation
	@param  double *vi - nodes of quadrature approximation
	@param  int &nNodes - number of nodes
	@param  int *mOrder - order of moments
	@param  double &CO2_l - weight fraction of CO2 in liquid
	@param  double &L_l - weight fraction of blowing agent in liquid
	@param  double &tmp - temperature
	@param  double &dVdt_2 - growth rate due to generation of CO2
	@param  double &dVdt_1 - growth rate due to evaporation of blowing agent
	@return void
*/

void growthSource(double *, double *, double *, double *, int &, int *, double &, double &, double &, double *, double *); // function prototype

void growthSource(double *sgBA, double *sgCO2, double *we, double *vi, int &nNodes, int *mOrder, double &CO2_l, double &L_l, double &tmp, double *volumeGrowthBA, double *volumeGrowthCO2)
{

	int i;
	int counter = 0;
	double k;

	while(counter<2*nNodes)
	{
		sgBA[counter]   = 0.0;
		sgCO2[counter]  = 0.0;

		// using static_cast to convert moment order from int to double
		k = static_cast<double>(mOrder[counter]);

		if(counter == 0)
		{
			sgBA[counter]  = 0.0;
			sgCO2[counter] = 0.0;
		}
		else if(counter == 1)
		{
			for(i=0;i<nNodes;i++)
			{
				if(vi[i] > 0.0)
				{
					sgBA[counter]  += volumeGrowthBA[i]*we[i];
					sgCO2[counter] += volumeGrowthCO2[i]*we[i];

				}
				else
				{
					sgBA[counter]   = sgBA[counter];
					sgCO2[counter]  = sgCO2[counter];
				}

			}

		}
		else
		{
			for(i=0;i<nNodes;i++)
			{
				if(vi[i] > 0.0)
				{
					sgBA[counter]  += k*volumeGrowthBA[i]*we[i]*(pow(vi[i], (k-1)));
					sgCO2[counter] += k*volumeGrowthCO2[i]*we[i]*(pow(vi[i], (k-1)));
				}
				else
				{
					sgBA[counter]   = sgBA[counter];
					sgCO2[counter]  = sgCO2[counter];
				}
			}
		}

		counter++;
	} // end of while
}
