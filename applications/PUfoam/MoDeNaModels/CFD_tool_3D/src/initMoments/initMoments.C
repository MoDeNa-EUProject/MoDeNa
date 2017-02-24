#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

extern "C" {void dsteqr_(char &, int *, double *, double *, double *, int *, double *, int *);}

#include "PDA.H"

int main(int argc, char *argv[])
{
	ifstream file("initMoments.txt");
	if(file.is_open())
	{
		int n;
		file >> n;
		n = n/2;
		double m[2*n];
		double we[n], xi[n];
		for (int i = 0; i < 2*n; i++)
		{
			file >> m[i];
		}
		PDA(we, xi, m, n);
		cout << "\n";

		for(int j=0;j<n;j++)
		{
			cout << "weight[" << j << "] = " << we[j] << endl;
		}
		for(int j=0;j<n;j++)
		{
			cout << "node[" << j << "] = " << xi[j] << endl;
		}

	}
	return 0;
}
