int HankelHadamard(const double *m, int nNodes);

int HankelHadamard(const double *m, int nNodes)
{
	int realizable = 0;
		
	if (nNodes == 2)
	{
		double** HH01 = new double*[nNodes];
		double** HH11 = new double*[nNodes];
		
		for(int i = 0; i < nNodes; ++i)
		{
			HH01[i] = new double[nNodes];
			HH11[i] = new double[nNodes];	
		}
		HH01[0][0]=m[0];
		HH01[0][1]=m[1];
		HH01[1][0]=m[1];
		HH01[1][1]=m[2];

		HH11[0][0]=m[1];
		HH11[0][1]=m[2];
		HH11[1][0]=m[2];
		HH11[1][1]=m[3];
		cout << endl;
		cout << "HH01 determinant: " << determinant(HH01, nNodes) << endl;
		cout << "HH11 determinant: " << determinant(HH11, nNodes) << endl;
		if (determinant(HH01, nNodes) < 0.0 || determinant(HH11, nNodes) < 0.0)
		{
			realizable = 1;
		}
	}
	
	if (nNodes == 3)
	{
		int twoNodes = nNodes - 1;
		double** HH01 = new double*[twoNodes];
		double** HH11 = new double*[twoNodes];
		for(int i = 0; i < twoNodes; ++i)
		{
			HH01[i] = new double[twoNodes];
			HH11[i] = new double[twoNodes];	
		}

		double** HH02 = new double*[nNodes];
		double** HH12 = new double*[nNodes];
		for(int i = 0; i < nNodes; ++i)
		{
			HH02[i] = new double[nNodes];
			HH12[i] = new double[nNodes];	
		}

		HH01[0][0]=m[0];
		HH01[0][1]=m[1];
		HH01[1][0]=m[1];
		HH01[1][1]=m[2];

		HH11[0][0]=m[1];
		HH11[0][1]=m[2];
		HH11[1][0]=m[2];
		HH11[1][1]=m[3];
		
		HH02[0][0]=m[0];
		HH02[0][1]=m[1];
		HH02[0][2]=m[2];
		HH02[1][0]=m[1];
		HH02[1][1]=m[2];
		HH02[1][2]=m[3];
		HH02[2][0]=m[2];
		HH02[2][1]=m[3];
		HH02[2][2]=m[4];
		
		HH12[0][0]=m[1];
		HH12[0][1]=m[2];
		HH12[0][2]=m[3];
		HH12[1][0]=m[2];
		HH12[1][1]=m[3];
		HH12[1][2]=m[4];
		HH12[2][0]=m[3];
		HH12[2][1]=m[4];
		HH12[2][2]=m[5];
		cout << endl;
		cout << "HH01 determinant: " << determinant(HH01, twoNodes) << endl;
		cout << "HH11 determinant: " << determinant(HH11, twoNodes) << endl;
		cout << "HH02 determinant: " << determinant(HH02, nNodes) << endl;
		cout << "HH12 determinant: " << determinant(HH12, nNodes) << endl;

		if (determinant(HH01, twoNodes) < 0.0 || determinant(HH11, twoNodes) < 0.0 || determinant(HH02, nNodes) < 0.0 || determinant(HH12, nNodes) < 0.0)
		{
			realizable = 1;
		}
	}

	if (nNodes == 4)
	{
		int twoNodes = nNodes - 2;
		int threeNodes = nNodes - 1;

		double** HH01 = new double*[twoNodes];
		double** HH11 = new double*[twoNodes];
		for(int i = 0; i < twoNodes; ++i)
		{
			HH01[i] = new double[twoNodes];
			HH11[i] = new double[twoNodes];	
		}

		double** HH02 = new double*[threeNodes];
		double** HH12 = new double*[threeNodes];
		for(int i = 0; i < threeNodes; ++i)
		{
			HH02[i] = new double[threeNodes];
			HH12[i] = new double[threeNodes];	
		}

		double** HH03 = new double*[nNodes];
		double** HH13 = new double*[nNodes];
		for(int i = 0; i < nNodes; ++i)
		{
			HH03[i] = new double[nNodes];
			HH13[i] = new double[nNodes];	
		}

		HH01[0][0]=m[0];
		HH01[0][1]=m[1];
		HH01[1][0]=m[1];
		HH01[1][1]=m[2];

		HH11[0][0]=m[1];
		HH11[0][1]=m[2];
		HH11[1][0]=m[2];
		HH11[1][1]=m[3];

		HH02[0][0]=m[0];
		HH02[0][1]=m[1];
		HH02[0][2]=m[2];
		HH02[1][0]=m[1];
		HH02[1][1]=m[2];
		HH02[1][2]=m[3];
		HH02[2][0]=m[2];
		HH02[2][1]=m[3];
		HH02[2][2]=m[4];

		HH12[0][0]=m[1];
		HH12[0][1]=m[2];
		HH12[0][2]=m[3];
		HH12[1][0]=m[2];
		HH12[1][1]=m[3];
		HH12[1][2]=m[4];
		HH12[2][0]=m[3];
		HH12[2][1]=m[4];
		HH12[2][2]=m[5];

		HH03[0][0]=m[0];
		HH03[0][1]=m[1];
		HH03[0][2]=m[2];
		HH03[0][3]=m[3];
		HH03[1][0]=m[1];
		HH03[1][1]=m[2];
		HH03[1][2]=m[3];
		HH03[1][3]=m[4];
		HH03[2][0]=m[2];
		HH03[2][1]=m[3];
		HH03[2][2]=m[4];
		HH03[2][3]=m[5];
		HH03[3][0]=m[3];
		HH03[3][1]=m[4];
		HH03[3][2]=m[5];
		HH03[3][3]=m[6];
		
		HH13[0][0]=m[1];
		HH13[0][1]=m[2];
		HH13[0][2]=m[3];
		HH13[0][3]=m[4];
		HH13[1][0]=m[2];
		HH13[1][1]=m[3];
		HH13[1][2]=m[4];
		HH13[1][3]=m[5];
		HH13[2][0]=m[3];
		HH13[2][1]=m[4];
		HH13[2][2]=m[5];
		HH13[2][3]=m[6];
		HH13[3][0]=m[4];
		HH13[3][1]=m[5];
		HH13[3][2]=m[6];
		HH13[3][3]=m[7];
		cout << endl;
		cout << "HH01 determinant: " << determinant(HH01, twoNodes) << endl;
		cout << "HH11 determinant: " << determinant(HH11, twoNodes) << endl;
		cout << "HH02 determinant: " << determinant(HH02, threeNodes) << endl;
		cout << "HH12 determinant: " << determinant(HH12, threeNodes) << endl;
		cout << "HH03 determinant: " << determinant(HH03, nNodes) << endl;
		cout << "HH13 determinant: " << determinant(HH13, nNodes) << endl;

		if ( determinant(HH01, twoNodes) < 0.0 || determinant(HH11, twoNodes) < 0.0 || determinant(HH02, threeNodes) < 0.0 || determinant(HH12, threeNodes) < 0.0 || determinant(HH03, nNodes) < 0.0 || determinant(HH13, nNodes) < 0.0 )
		{
			realizable = 1;
		}
	}
	return realizable;
}