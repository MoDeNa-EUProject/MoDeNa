void enterMoments(double *mom, const int nNodes);
void enterWeightsNodes(double *we, double *vi, int nNodes);
void weightsNodesToMoms(double *mom, double *we, double *vi, int nNodes);
void printMoms(const double *m, int nNodes);
bool momentsPositivite(const double *moments, int &nNodes);
void normalizeMom(double *m, int nNodes);
void denormalizeMom(double *m, double M0, int nNodes);
void initializeDifferenceTable(double **differenceTable, int nNodes);
void printDifferenceTable(double **differenceTable, int nNodes);
void buildDifferenceTable(double **differenceTable, double *mom, int nNodes);
int isRealizable(double **differenceTable, int nNodes);
double vectorNorm(vector<double> const& v);

void enterMoments(double *mom, const int nNodes)
{
	for (int i = 0; i < 2*nNodes; i++)
	{
		cout << "\tm[" << i << "] = "; cin >> mom[i];		
	}
	cout << endl;
}
void enterWeightsNodes(double *we, double *vi, int nNodes)
{
	for(int j = 0; j < nNodes; j++)
		{
			cout << "\twe[" << j << "] = "; cin >> we[j];
			cout << "\tvi[" << j << "] = "; cin >> vi[j];
			cout << endl;
		}
}
void weightsNodesToMoms(double *mom, double *we, double *vi, int nNodes)
{
	for(int k = 0; k < 2*nNodes; k++)
	{
		for (int i = 0; i < nNodes; i++)
		{
			if (k == 0)
			{
				mom[k] += we[i];
			}
			else
			{
				mom[k] += we[i]*pow(vi[i],k);			
			}
		}
	}
}
void printMoms(const double *m, int nNodes)
{
	for(int i = 0; i < 2*nNodes; i++)
	{
		cout << "\tm[" << i << "] = " << m[i] << endl;
	}
}

bool momentsPositivite(const double *moments, int &nNodes)
{
	int noMoms = 2*nNodes;
	for (int i = 0; i < noMoms; i++)
	{
		if (moments[i] < 0.0)
		{
			cout << "m[" << i << "] is negative!" << endl;
			return false;
			// exit(1);
		}	
	}
	return true;
}

void normalizeMom(double *m, int nNodes)
{
	for (int i = 2*nNodes; i >= 0; i--)
	{
		m[i] = m[i]/m[0];
	}
}

void denormalizeMom(double *m, double M0, int nNodes)
{
	for (int i = 0; i < 2*nNodes; i++)
	{
		m[i] = m[i]*M0;
	}
}

void initializeDifferenceTable(double **differenceTable, int nNodes)
{
	for (int r = 0; r < 2*nNodes; r++)
    {
    	for (int c = 0; c < 2*nNodes; c++)
    	{
    		differenceTable[r][c] = 0.0;
    	}
    }
}

void printDifferenceTable(double **differenceTable, int nNodes)
{
	cout << endl << "Table of Difference: \n";
	for (int r = 0; r < 2*nNodes; r++)
    {
    	for (int c = 0; c < 2*nNodes; c++)
    	{
    		// formatting decimals
    		cout.setf(ios::fixed);
			cout.setf(ios::showpoint);
			cout.precision(10);
    		cout << differenceTable[r][c] << " ";
    	}
    	cout << endl;
    } 
}

void buildDifferenceTable(double **differenceTable, double *mom, int nNodes)
{
	for (int i = 0; i < 2*nNodes; i++)
	{
		differenceTable[i][0] = log(mom[i]);
	}
	for (int j = 1; j <= 2*nNodes; j++) 
    {
        for (int k = 0; k < (2*nNodes - j); k++)
        {   
        	differenceTable[k][j] = differenceTable[k+1][j-1]-differenceTable[k][j-1];    
        }            
    }
}

int isRealizable(double **differenceTable, int nNodes)
{
	int realizable = 0;
	for (int r = 0; r < (2*nNodes-2); r++)
	{
		if (differenceTable[r][2] < 0.0)
		{
			// cerr << "Warning Unrealizable Moments!";
			realizable = 1;
		}
	}
	return realizable;
}

double vectorNorm(vector<double> const& v)
{
	double result = 0.0;
	for (unsigned int i = 0; i < v.size(); ++i)
	{
		result += v[i]*v[i];
	}
	return (sqrt(result));
}
