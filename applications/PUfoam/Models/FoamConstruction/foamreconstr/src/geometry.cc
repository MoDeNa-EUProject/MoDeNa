/*! \file
	\brief Helpful geometric and mathematical functions.
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#include "globals.hh"
using namespace globals;
//! True if the point is in the domain, False otherwise.
bool in_domain(\
	double x /**< [in] `x` position */, \
	double y /**< [in] `x` position */, \
	double z /**< [in] `x` position */)
{
	if (x>0 && x<nx && y>0 && y<ny && z>0 && z<nz)
		return true;
	else
		return false;
}
//! Cross product of two points in 3D.
double *CrossProduct(\
	double a[] /**< [in] 1st point */,\
	double b[] /**< [in] 2nd point */)
{
    double *cp;
    cp = new double[3];
    cp[0]=a[1]*b[2]-a[2]*b[1];
    cp[1]=a[2]*b[0]-a[0]*b[2];
    cp[2]=a[0]*b[1]-a[1]*b[0];
    return(cp);
}
//! True if two points are in the same half-plane defined by other two points.
bool SameSide(\
	double p1[] /**< [in] 1st point */,\
	double p2[] /**< [in] 2nd point */,\
	double a[] /**< [in] 1st point defining half-plane */,\
	double b[] /**< [in] 2nd point defining half-plane */)
{
    bool dec;
    double ba[3],p1a[3],p2a[3];
    int i;

    for (i=0; i<3; i++) {
        ba[i]=b[i]-a[i];
        p1a[i]=p1[i]-a[i];
        p2a[i]=p2[i]-a[i];
    }
    double *cp1=CrossProduct(ba,p1a);
    double *cp2=CrossProduct(ba,p2a);
    dec=false;
    if ((cp1[0]*cp2[0]+cp1[1]*cp2[1]+cp1[2]*cp2[2])>=0) dec=true;
    delete cp1;
    delete cp2;
    return(dec);
}
//! True if point lies in the triangle.
bool PointInTriangle(\
	double p[] /**< [in] the point */,\
	double a[] /**< [in] 1st point defining triangle */,\
	double b[] /**< [in] 2nd point defining triangle */,\
	double c[] /**< [in] 3rd point defining triangle */)
{
    bool dec;
    dec=false;
    if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b)) {
		dec=true;
	}
    return(dec);
}
//! Calculates porosity.
double porosity(int ***amat /**< [in] matrix */) {
	double por=0;
	int i,j,k;
	for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                if (amat[i][j][k]==0) por++;
            }
        }
    }
	por /= nx*ny*nz;
	return(por);
}
