/**
@ingroup mod_0Dcfd
@file initializeMoments.h
@brief initializes the moments
@fn void mom_init(double *momz, double &size_0, int &n, double &sigma, double &NN)
@param momz moments of different orders
@param size_0 initial bubble diameter
@param n number of moments
@param sigma correlated to the variance of initial distribution
@param NN correlated to the number of initial bubbles per unit volume
@return void
*/
void mom_init(double *, double &, int &, double &, double &);

void mom_init(double *momz, double &size_0, int &n, double &sigma, double &NN)
{
	int i;
	double v_0 	= (M_PI/6.0)*(pow(size_0,3));
	double mu 	= log(v_0);

	for(i=0;i<n;i++)
	{
		if(i == 0)
		{
			momz[0] = NN;
		}
		else
		{
			momz[i] = NN*exp(i*mu+0.5*i*i*sigma*sigma);
		}

	}
}