/**
@ingroup mod_0Dcfd
@file differenceTable.h
@brief Calculate difference table for the moments realizability check
@fn void enterMoments(double *mom, const int nNodes)
@brief input the moments
@param mom moments of different orders
@param nNodes number of nodes
@return void
@fn void enterWeightsNodes(double *we, double *vi, int nNodes)
@brief inputs weights and nodes of quadrature
@param we weights of quadrature
@param vi nodes of quadrature
@param nNodes number of nodes
@return void
@fn void weightsNodesToMoms(double *mom, double *we, double *vi, int nNodes);
@brief converts weights and nodes to moments
@param mom moments of different orders
@param we weights of quadrature
@param vi nodes of quadrature
@param nNodes number of nodes
@return void
@fn void printMoms(const double *m, int nNodes)
@brief prints moments
@param m moments of different orders
@param nNodes number of nodes
@return void
@fn bool momentsPositivite(const double *moments, int &nNodes)
@param moments moments of different orders
@param nNodes number of nodes
@return true if moments are positive
@fn void normalizeMom(double *m, int nNodes)
@brief normalizes the moments
@param m moments of different orders
@param nNodes number of nodes
@return void
@fn void denormalizeMom(double *m, double M0, int nNodes)
@brief de-normalizes the moments
@param m moments of different orders
@param M0 moment of order zero based on unit volume of the foam
@param nNodes number of nodes
@return void
@fn void initializeDifferenceTable(double **differenceTable, int nNodes)
@brief initializes the difference table
@param differenceTable pointer to difference table
@param nNodes number of nodes
@return void
@fn void printDifferenceTable(double **differenceTable, int nNodes)
@brief prints the difference table
@param differenceTable pointer to difference table
@param nNodes number of nodes
@return void
@fn void buildDifferenceTable(double **differenceTable, double *mom, int nNodes)
@brief constructs the difference table
@param differenceTable pointer to difference table
@param mom moments of different orders
@param nNodes number of nodes
@return void
@fn int isRealizable(double **differenceTable, int nNodes)
@brief checks if the moments are realizable
@param differenceTable pointer to the difference table
@param nNodes number of nodes
@return 1 if the moments are realizable
@fn double vectorNorm(vector<double> const& v)
@brief normalized vector used in McGraw correction algorithm
@param v input vector
@return normalized vector
*/
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
