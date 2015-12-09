/** @file coalescence.h 
	@brief Source term due to bubble coalescence 
	@fn void coalescenceSource(double *sc, double *we, double *vi, int &nNodes, int *mOrder, double &beta0)
	@param  double *sc - source due to coalescence
	@param  double *we - weights of quadrature approximation
	@param  double *vi - nodes of quadrature approximation
 	@param  int &nNodes - number of nodes
	@param  int *mOrder - order of moments
	@param  double &beta0 - coalescence constant
	@return void. 
*/
	
void coalescenceSource(double *, double *, double *, int &, int *, double &);

void coalescenceSource(double *sc, double *we, double *vi, int &nNodes, int *mOrder, double &beta0)
{

	
	int counter = 0;	
	double k; 			// To hold the moment order
	int i, j;
	
// coalescence kernel beta0*(v(i) + v(j))

	while(counter<2*nNodes)
	{	
		sc[counter] = 0.0;
		k = static_cast<double>(mOrder[counter]);

		if(counter == 0)
		{	
			for(i=0;i<nNodes;i++)
			{
				for(j=0;j<nNodes;j++)
				{	
					if (vi[i]*vi[j] != 0)
					{
						sc[counter] += 0.5*beta0*(vi[i]+vi[j])*we[i]*we[j]*(-1.0);
					}
					else
					{
						sc[counter] = 0.0;
					}	
				}
			}			
		}
		else if(counter == 1)
		{
			sc[counter] = 0.0;
		}		
		else
		{
			for(i=0;i<nNodes;i++)
			{
				for(j=0;j<nNodes;j++)
				{	
					if(vi[i]*vi[j] != 0)
					{
						sc[counter] += 0.5*beta0*(vi[i]+vi[j])*we[i]*we[j]*(pow((vi[i]+vi[j]),k)-pow(vi[i],k)-pow(vi[j],k));
					}
					else
					{
						sc[counter]  = 0.0;	
					}
				}
			}	
		}
		counter++;
	} // End of while loop
			
}

/*Binomial coefficient or kChoosen 
int choose(int, int);

int choose(int k, int n)
{
	if(n == 0)
	{
		return 1;
	}
	else
	{
		return (k*choose(k-1, n-1))/k;
	} 

}*/
