void WrightCorrection(double *m, int nNodes);

void WrightCorrection(double *m, int nNodes)
{
	int n = 2*nNodes;
	// dynamic allocation of new moments
	double *m1, *m2;
	m1 = new double[n];
	m2 = new double[n]; 

	// selecting two moments, 2nd and 3rd moments
	int i = 2; 
	int j = 3; 
	
	// corresponding mean and variance
	double mu =  (j/(i*j-pow(i,2)))*log(m[i]/m[0]) + (i/(i*j-pow(j,2)))*log(m[j]/m[0]);
	double sigma_var = ((2.0/(pow(j,2)))*log(m[j]/m[0])-(2.0/(i*j))*log(m[i]/m[0]))/(1.0-i/j);
	if (sigma_var < 0) 
	{	
		sigma_var = 0.0;
	}

	// new moments for the first distribution
	for (int k = 0; k < n; k++)
	{
		m1[k] = m[0]*exp(k*mu+((pow(k,2)*sigma_var)/2.0));
	}

	// selecting two more moments
	i = 1;
	j = 3;

	// corresponding mean and variance
	mu =  (j/(i*j-pow(i,2)))*log(m[i]/m[0]) + (i/(i*j-pow(j,2)))*log(m[j]/m[0]);
	sigma_var = ((2.0/(pow(j,2)))*log(m[j]/m[0])-(2.0/(i*j))*log(m[i]/m[0]))/(1.0-i/j);
	if (sigma_var<0) 
	{
    	sigma_var = 0.0;
	}

	// new moments for the second distribution
	for (int k = 0; k < n; k++)
	{
		m2[k] = m[0]*exp(k*mu+((pow(k,2)*sigma_var)/2.0));
	}

	// corrected moments = arithmetic mean of the previous two moment sets.
	for (int k = 0; k < n; k++)
	{
		m[k] = (m1[k] + m2[k])/2.0;
	}
}