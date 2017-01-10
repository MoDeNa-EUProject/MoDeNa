/**
@ingroup mod_3Dcfd
@file coalescenceSource.H
@brief source term due to bubbles coalescence
@fn void coalescenceSource(double *sc, double *we, double *vi, int &nNodes, int *mOrder)
@param sc source of coalescence
@param we weights of quadrature
@param vi nodes of quadrature
@param nNodes number of nodes
@param mOrder order of moments
@returns void
*/
void coalescenceSource(double *, double *, double *, int &, int *);

void coalescenceSource
(
    double *sc, double *we, double *vi, int &nNodes, int *mOrder
)
{
    int counter = 0;
    double k;
    int i, j;
    double coalescenceRate;
    while (counter < 2*nNodes)
    {
        sc[counter] = 0.0;
        k = static_cast<double>(mOrder[counter]);

        coalescenceRate = coalescenceKernel();

        if (counter == 0)
        {
            for (i = 0; i < nNodes; i++)
            {
                for (j = 0; j < nNodes; j++)
                {
                    if (vi[i]*vi[j] != 0)
                    {
                        sc[counter] +=
                        (
                            coalescenceRate*(vi[i] + vi[j])*we[i]*we[j]*(-1.0)
                        );
                    }
                    else
                    {
                        sc[counter] = 0.0;
                    }
                }
            }
        }
        else if (counter == 1)
        {
            sc[counter] = 0.0;
        }
        else
        {
            for (i = 0; i < nNodes; i++)
            {
                for (j = 0; j < nNodes; j++)
                {
                    if (vi[i]*vi[j] != 0)
                    {
                        sc[counter] +=
                        (
                            coalescenceRate*(vi[i] + vi[j])*we[i]*we[j]
                            *(Foam::pow((vi[i] + vi[j]),k)
                            - Foam::pow(vi[i],k) - Foam::pow(vi[j],k))
                        );
                    }
                    else
                    {
                        sc[counter]  = 0.0;
                    }
                }
            }
        }
        counter++;
    }
}