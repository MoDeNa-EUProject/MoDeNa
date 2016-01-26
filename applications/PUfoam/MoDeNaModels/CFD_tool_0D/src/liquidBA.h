/** @file liquidBA.h
	@brief Maximum soluble physical blowing agent in the liquid mixture
	@fn double LMax (double tm)
	@param double tm - the input temperature 
	@return maximum soluble blowing agent
*/
double LMax(double);

double LMax (double tm)
{
	double lMax,aa,hh,T0,ww;
	switch (phBL) {
		case 1:
			aa 		= 0.0064;	// -
			hh 		= 0.0551;	// -
			T0 		= 298.0;	// K
			ww 		= 17.8;		// 1/K
			// double T_in 	= 300; 		// K
			if (tm > T0)
			{
				double tempDummy = pow((tm-T0),2.0);
				lMax = (aa + hh*exp((-tempDummy/(2.0*ww*ww))));
				break;
			}
			else
			{
				lMax=aa+hh;
				break;
			}
		case 2:
			aa 		= 1e-7;	// -
			hh 		= 4.2934;	// -
			T0 		= 203.3556;	// K
			ww 		= 40.016;		// 1/K
			// double T_in 	= 300; 		// K
			if (tm > T0)
			{
				double tempDummy = pow((tm-T0),2.0);
				lMax = (aa + hh*exp((-tempDummy/(2.0*ww*ww))));
				break;
			}
			else
			{
				lMax=aa+hh;
				break;
			}
	}
	return lMax;
}
