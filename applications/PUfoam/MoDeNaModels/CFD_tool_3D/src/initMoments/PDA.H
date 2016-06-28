// inline bool KroneckerDelta(const int i, const int j)
// {
// return (i==j);
// }


void PDA(double *, double *, double *, int &);


void PDA(double *we, double *vi, double *mom, int &n)
{

// Construct P Matrix 
    int i,j;

    double p[2*n+1][2*n+1];

    for(i=0;i<2*n+1;i++)
    {
        for(j=0;j<2*n+1;j++)
        {
            p[i][j] = 0.0;
        }
    }   
    
    double norm_mom[2*n];
   
// Normalize the moments
    for(i=0;i<2*n;i++)
    {
        // if(mom[0] < 1.0e-20)
        // {
        //     norm_mom[i] = mom[i];
        // }
        // else
        // {
            norm_mom[i] = mom[i]/(max(mom[0],1.0e-10));
        // }
    }	
// First column of P matrix
    // for(i=0;i<2*n;i++)
    // {
    //     p[i][0] = KroneckerDelta(i+1,1);
    // }
    p[0][0] = 1.0;
    p[0][1] = 1.0;
// Second column of P matrix
    for(i=1;i<2*n;i++)
    {
        p[i][1] = pow(-1,double (i))*norm_mom[i]; // Foam::pow 
    }

// Recursion method for calculation of P
    for(j=2;j<2*n+1;j++)
    {
        for(i=0;i<2*n+2-j;i++)
            {
                p[i][j] = p[0][j-1]*p[i+1][j-2] - p[0][j-2]*p[i+1][j-1];
            }
    }    

// Computing zeta 
    double zeta[2*n];
    for(i=0;i<2*n;i++)
    {
        zeta[i] = 0.0;
    }

    for(i=1;i<2*n;i++)
    {
        if(p[0][i]*p[0][i-1]>0)
        {
            zeta[i] = p[0][i+1]/(p[0][i]*p[0][i-1]);
        }
	    else
	    {  
	        zeta[i] = 0.0;
	    }
    }

// Coefficients of Jacobi matrix
    double aa[n];
    double bb[n-1], cc[n-1];

    for(i=0;i<n-1;i++)
    {
        bb[i] = cc[i] = 0.0;
    }

    for(i=1;i<=n;i++)
    {
        aa[i-1]=zeta[2*i-1]+zeta[2*i-2];
    }
    
    for(i=1;i<=n-1;i++)
    {
        bb[i-1]=zeta[2*i]*zeta[2*i-1];
    }
    
    for(i=1;i<=n-1;i++)
    {
        cc[i-1]= sqrt(fabs(bb[i-1]));
    }


// Eigenvalues and Eigenvectors of Jacobi matrix
    double evec[n][n];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            evec[i][j] = 0.0;
        }
    }


    double work[2*n-2];
    int info;
    char choice='I';
	
    dsteqr_(choice,&n,aa,cc,&evec[0][0],&n,work,&info);
	
    for(i=0;i<n;i++)
    {
        vi[i]=aa[i];
        
        if(vi[i] < 0.0)
        {
            vi[i] = 0.0;
        }
	}

    for(j=0;j<n;j++)
    {
        we[j]=mom[0]*pow(evec[j][0],2); // Foam::pow
        
        if(we[j] < 0)
        {
            we[j] = 0.0;
        }

        // cout<< "we[" << j << "] = " << we[j] << ", vi[" << j << "] = " << vi[j] << endl;
    }
} // End of PDA

