/*! \file
	\brief Functions for creating struts at cell vertices.
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#include "globals.hh"
#include <iostream>
#include <cmath>
#include "determinant.hh"
#include "geometry.hh"
#include "allocation.hh"
using namespace std;
using namespace globals;
//! Function for creating struts at cell vertices.
//!
//! Struts are in shape of tetrahedra. It is assumed that four cells
//! share a vertex. Vertices are taken from tessellation. Size of the
//! tetrahedron is determined by the `dstrut` parameter. Changes `amat`.
void makeNodeStruts(\
    int ***amat /**< [in,out] voxel matrix */,\
    int sv /**< [in] number of vertices */,\
    int incmax /**< [in] maximum number of vertex connections */,\
    int vmax /**< [in] maximum number of vertices */,\
    double **vert /**< [in] vertex positions */,\
    int **vinc /**< [in] indexes of connected vertices */,\
    bool report /**< [in] show output */)
{
    // Creates strut parts - solid material around cell vertices. Created
    // strut has a shape of tetrahedron. Updates `amat`.
    int i,j,k,l,m,n;
    int svd; //final number of vertices in the domain
    double xx;
    double mini[3],maxi[3];
    double xmin=0,xmax=nx;
    double ymin=0,ymax=ny;
    double zmin=0,zmax=nz;
    double dist[4];
    double ***tetra; //four points defining each tetrahedron
	tetra=alloc_3Ddmatrix(vmax,4,3);
    // calculate tetrahedra
    if (report) {
        cout << "creating struts at cell vertices..." << endl;
    }
    svd=0;
    for (i=0; i<sv; i++) {
        if (in_domain(vert[i][0],vert[i][1],vert[i][2])) {
            k=0;
            for (j=0; j<incmax; j++) {
                if (vinc[i][j] != -1) {
                    k++;
                } else {
                    break;
                }
            }
            if (k>3) {
                for (l=0; l<4; l++) {
                    xx=0;
                    for (m=0; m<3; m++) {
                        xx += pow(vert[vinc[i][l]][m]-vert[i][m],2);
                    }
                    dist[l]=sqrt(xx);
                }
                for (l=0; l<4; l++) {
                    for (m=0; m<3; m++) {
                        if (dstrut/dist[l]<1) {
                            tetra[svd][l][m]=vert[i][m] + \
                                dstrut/dist[l]*(vert[vinc[i][l]][m]-\
                                vert[i][m]); //strut parameter
                        } else {
                            tetra[svd][l][m]=vert[vinc[i][l]][m];
                            //for safety
                            //(we don't want sharp edges (just try and see))
                        }
                    }
                }
                svd++;
            }
        }
    }
    // determine if voxel is inside the tetrahedron
    double **A0,**A1,**A2,**A3,**A4;
    A0 = new double*[4]; //allocate
    for (i=0; i<4; ++i)
    A0[i] = new double[4];
    A1 = new double*[4]; //allocate
    for (i=0; i<4; ++i)
    A1[i] = new double[4];
    A2 = new double*[4]; //allocate
    for (i=0; i<4; ++i)
    A2[i] = new double[4];
    A3 = new double*[4]; //allocate
    for (i=0; i<4; ++i)
    A3[i] = new double[4];
    A4 = new double*[4]; //allocate
    for (i=0; i<4; ++i)
    A4[i] = new double[4];
    for (i=0; i<svd; i++) {
        mini[0]=2*xmax; maxi[0]=-2*xmax;
        mini[1]=2*ymax; maxi[1]=-2*ymax;
        mini[2]=2*zmax; maxi[2]=-2*zmax;
        for (j=0; j<4; j++) {
            for (k=0; k<3; k++) {
                if (tetra[i][j][k]<mini[k]) mini[k]=tetra[i][j][k];
                if (tetra[i][j][k]>maxi[k]) maxi[k]=tetra[i][j][k];
            }
        }
        for (j=(int) floor(mini[0]); j<=ceil(maxi[0]); j++)
            for (k=(int) floor(mini[1]); k<=ceil(maxi[1]); k++)
                for (l=(int) floor(mini[2]); l<=ceil(maxi[2]); l++) {
                    for (m=0; m<4; m++) {
                        for (n=0; n<3; n++) {
                            A0[m][n]=tetra[i][m][n];
                            A1[m][n]=tetra[i][m][n];
                            A2[m][n]=tetra[i][m][n];
                            A3[m][n]=tetra[i][m][n];
                            A4[m][n]=tetra[i][m][n];
                        }
                    }
                    for (m=0; m<4; m++) {
                        A0[m][3]=1;
                        A1[m][3]=1;
                        A2[m][3]=1;
                        A3[m][3]=1;
                        A4[m][3]=1;
                    }
                    A1[0][0]=j; A1[0][1]=k; A1[0][2]=l;
                    A2[1][0]=j; A2[1][1]=k; A2[1][2]=l;
                    A3[2][0]=j; A3[2][1]=k; A3[2][2]=l;
                    A4[3][0]=j; A4[3][1]=k; A4[3][2]=l;
                    m=sgn(Determinant(A0,4));
                    if (m==sgn(Determinant(A1,4)) && \
                        m==sgn(Determinant(A2,4)) && \
                        m==sgn(Determinant(A3,4)) && \
                        m==sgn(Determinant(A4,4))) {
                        amat[CONFINEX(j)][CONFINEY(k)][CONFINEZ(l)]=1;
                    }
                }
    }
    tetra=free_3Ddmatrix(tetra);
    for (i=0; i<4; ++i) delete [] A0[i]; //deallocate
    delete [] A0;
    for (i=0; i<4; ++i) delete [] A1[i]; //deallocate
    delete [] A1;
    for (i=0; i<4; ++i) delete [] A2[i]; //deallocate
    delete [] A2;
    for (i=0; i<4; ++i) delete [] A3[i]; //deallocate
    delete [] A3;
    for (i=0; i<4; ++i) delete [] A4[i]; //deallocate
    delete [] A4;
    return;
}
