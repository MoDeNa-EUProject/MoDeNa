// source term due to bubbles coalescence
// @param - double *sc - source of coalescence
// @param - double *we - weights of quadrature
// @param - duoble *vi - nodes of quadrature
// @param - int &nNodes - number of nodes
// @param - int *mOrder - order of moments
// Returns - void
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