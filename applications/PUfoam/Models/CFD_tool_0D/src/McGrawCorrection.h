/**
@file McGrawCorrection.h
@brief McGraw correction algorithm 
@fn double cos_quad_alpha(double **bk, double **differenceTable, int k, int nNodes);
@brief 
*/
double cos_quad_alpha(double **bk, double **differenceTable, int k, int nNodes);
void McGrawCorrection(double *moments, int nNodes);


double cos_quad_alpha(double **bk, double **differenceTable, int k, int nNodes)
{
    int n = 2*nNodes;
    std::vector<double> bk_kthRow(n);
    std::vector<double> dt_3rdColumn(n);
    for (int i = 0; i < n; i++)
    {
        bk_kthRow[i]    = bk[k][i];
        dt_3rdColumn[i] = differenceTable[i][3]; 
    }
    double nominator   = inner_product(bk_kthRow.begin(), bk_kthRow.end(), dt_3rdColumn.begin(), 0.0);
    double denominator = vectorNorm(bk_kthRow)*vectorNorm(dt_3rdColumn);
    return ( pow((nominator/denominator),2) );
}

void McGrawCorrection(double *moments, int nNodes)
{
	int n = 2*nNodes;
	int maxIteration = 100;

	// dynamic allocation of bkk, bk and difference teble
	double** bkk 			 = new double*[n];
	double** bk  			 = new double*[n];
	double** differenceTable = new double*[n];
	
	for(int i = 0; i < n; ++i)
	{
		bkk[i] 				= new double[n];
		bk[i]  				= new double[n];
		differenceTable[i] 	= new double[n];	
	}

	// build bkk and bk
	for (int r = 0; r < n; r++)
	{
		for (int c = 0; c < n; c++)
		{
			bkk[r][c] = 0.0;
		}
	}
	for (int k = 0; k < n; k++)
    {   
        for (int i = 0; i < n; i++)
        {
            if (i == k)
            {
            	bkk[i][0] = 1.0 ;
            }
            else
            {
                bkk[i][0] = 0.0;
            }
        }
        
        for (int j=1; j <= 3; j++)
        {
            for (int i = 0; i < (n-j); i++)
            {
                bkk[i][j]=bkk[i+1][j-1]-bkk[i][j-1];
            }
        }
        // build bk        
        for (int i = 0; i < n; i+=1)
        {
            bk[k][i] = bkk[i][3];
        }
    }

    int iteration 	= 0;
    int realizable 	= 1;

    while (realizable == 1 && iteration < maxIteration)
    {
    	realizable = 0;
    	// initializeDifferenceTable(differenceTable, nNodes);
    	buildDifferenceTable(differenceTable, moments, nNodes);
    	realizable = isRealizable(differenceTable, nNodes);

    	if ( realizable == 1 )
    	{
    		iteration++;
    		// cout << "iteration: " << iteration << endl;
    		int k_star = 1;
    		
    		for (int k = 0; k < n; k++)
    		{
    			double cqa 		= cos_quad_alpha(bk, differenceTable, k, nNodes);
    			double cqa_star = cos_quad_alpha(bk, differenceTable, k_star, nNodes);
    			if(cqa >= cqa_star)
    			{
    				k_star = k;
    			}   			
    		}
    		std::vector<double> bk_starRow(n);
    		std::vector<double> dt_3rdColumn(n);
    		for (int i = 0; i < n; i++)
    		{
    			bk_starRow[i] 	= bk[k_star][i];
    			dt_3rdColumn[i] = differenceTable[i][3];
    		}
    		double lnck = -( inner_product(bk_starRow.begin(), bk_starRow.end(), dt_3rdColumn.begin(), 0.0) )/( vectorNorm(bk_starRow) );
    		moments[k_star] = exp(lnck)*moments[k_star];   		
    	}
    }
}