#include "globals.hh"
using namespace globals;
bool in_domain(float x, float y, float z) {
	if (x>0 && x<nx && y>0 && y<ny && z>0 && z<nz)
		return true;
	else
		return false;
}

float *CrossProduct(float a[],float b[]) {
    float *cp;
    cp = new float[3];
    cp[0]=a[1]*b[2]-a[2]*b[1];
    cp[1]=a[2]*b[0]-a[0]*b[2];
    cp[2]=a[0]*b[1]-a[1]*b[0];
    return(cp);
}

bool SameSide(float p1[],float p2[],float a[],float b[]) {
    bool dec;
    float ba[3],p1a[3],p2a[3];
    int i;

    for (i=0; i<3; i++) {
        ba[i]=b[i]-a[i];
        p1a[i]=p1[i]-a[i];
        p2a[i]=p2[i]-a[i];
    }
    float *cp1=CrossProduct(ba,p1a);
    float *cp2=CrossProduct(ba,p2a);
    dec=false;
    if ((cp1[0]*cp2[0]+cp1[1]*cp2[1]+cp1[2]*cp2[2])>=0) dec=true;
    delete cp1;
    delete cp2;
    return(dec);
}

bool PointInTriangle(float p[],float a[],float b[],float c[]) {
    bool dec;
    dec=false;
    if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b)) {
		dec=true;
	}
    return(dec);
}

double porosity(int ***amat) {
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
