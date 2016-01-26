/*! \file
	\brief Functions for creating struts at cell edges.
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#include "globals.hh"
#include <iostream>
#include <cmath>
#include "geometry.hh"
using namespace std;
using namespace globals;
//! Function for creating struts at cell edges.
//!
//! Struts are in shape of triangular prisms. It is assumed that three cells
//! share an edge. Edges are taken from tessellation. Size of the triangular
//! base is determined by the `dedge` parameter. Changes `amat`.
void makeEdgeStruts(\
    int ***amat /**< [in,out] voxel matrix */, \
    int sv /**< [in] number of vertices */, \
    double **vert /**< [in] vertex positions */, \
    int **vinc /**< [in] indexes of connected vertices */, \
    bool report /**< [in] show output */)
{
    // Creates struts - solid material along cell edge. Created strut has
    // a shape of triangular prism. Updates `amat`.
    int i,j,k,l,m,n;
    double da,db,dc,de; //distance
    double xx,xxa,xxb,xxc;
    double mini[3],maxi[3];
    double xmin=0,xmax=nx;
    double ymin=0,ymax=ny;
    double zmin=0,zmax=nz;
    double n0[3]; //primary edge unit vector
    double a0[3],b0[3],c0[3]; //incident edges unit vectors
    //(their projection to plane normal to primary edge vector)
    double pa[3],pb[3],pc[3]; //vertices of triangle base of the prism
    double pap[3],pbp[3],pcp[3]; //vertices of triangle in plane of point P
	double pp[3]; //coordinates of point P
    if (report) {
        cout << "creating struts at cell edges..." << dedge << endl;
    }
    for (i=0; i<sv; i++) { // for all vertices
        if (in_domain(vert[i][0],vert[i][1],vert[i][2])) {
            for (j=0; j<4; j++) { // for every edge
                de=sqrt(pow(vert[i][0]-vert[vinc[i][j]][0],2)+\
                    pow(vert[i][1]-vert[vinc[i][j]][1],2)+\
                    pow(vert[i][2]-vert[vinc[i][j]][2],2));
                    //length of the primary edge
                for (k=0; k<3; k++) {
                    // define primary and incident edge vectors
                    n0[k]=vert[vinc[i][j]][k]-vert[i][k];
                    a0[k]=vert[vinc[i][(j+1)%4]][k];
                    b0[k]=vert[vinc[i][(j+2)%4]][k];
                    c0[k]=vert[vinc[i][(j+3)%4]][k];
                }
                xx=0; // normalize primary vector
                for (k=0; k<3; k++) {
                    xx+=pow(n0[k],2);
                }
                for (k=0; k<3; k++) {
                    n0[k]/=sqrt(xx);
                }
                da=abs(n0[0]*a0[0]+n0[1]*a0[1]+n0[2]*a0[2]-\
                    n0[0]*vert[i][0]-n0[1]*vert[i][1]-n0[2]*vert[i][2]);
                    // calculate distance from plane normal to
                    // primary vector at cell vertex
                for (k=0; k<3; k++) {
                    // redefine incident edge vectors to be in plane normal
                    // to primary edge vector (project them to the plane)
                    a0[k]+=da*n0[k];
                }
                if (abs(n0[0]*a0[0]+n0[1]*a0[1]+n0[2]*a0[2]-\
                    n0[0]*vert[i][0]-n0[1]*vert[i][1]-\
                    n0[2]*vert[i][2])>1e-2) {
                    for (k=0; k<3; k++) {
                        a0[k]-=2*da*n0[k];
                    }
                }
                db=abs(n0[0]*b0[0]+n0[1]*b0[1]+n0[2]*b0[2]-\
                    n0[0]*vert[i][0]-n0[1]*vert[i][1]-n0[2]*vert[i][2]);
                for (k=0; k<3; k++) {
                    b0[k]+=db*n0[k];
                }
                if (abs(n0[0]*b0[0]+n0[1]*b0[1]+n0[2]*b0[2]-\
                    n0[0]*vert[i][0]-n0[1]*vert[i][1]-\
                    n0[2]*vert[i][2])>1e-2) {
                    for (k=0; k<3; k++) {
                        b0[k]-=2*db*n0[k];
                    }
                }
                dc=abs(n0[0]*c0[0]+n0[1]*c0[1]+n0[2]*c0[2]-\
                    n0[0]*vert[i][0]-n0[1]*vert[i][1]-n0[2]*vert[i][2]);
                for (k=0; k<3; k++) {
                    c0[k]+=dc*n0[k];
                }
                if (abs(n0[0]*c0[0]+n0[1]*c0[1]+n0[2]*c0[2]-\
                    n0[0]*vert[i][0]-n0[1]*vert[i][1]-\
                    n0[2]*vert[i][2])>1e-2) {
                    for (k=0; k<3; k++) {
                        c0[k]-=2*dc*n0[k];
                    }
                }
                for (k=0; k<3; k++) {
                    // redefine incident edge vectors to be in plane normal
                    // to primary edge vector (project them to the plane)
                    a0[k]-=vert[i][k];
                    b0[k]-=vert[i][k];
                    c0[k]-=vert[i][k];
                }
                xxa=0; xxb=0; xxc=0; // normalize incident edge vectors
                for (k=0; k<3; k++) {
                    xxa+=pow(a0[k],2);
                    xxb+=pow(b0[k],2);
                    xxc+=pow(c0[k],2);
                }
                for (k=0; k<3; k++) {
                    a0[k]/=sqrt(xxa);
                    b0[k]/=sqrt(xxb);
                    c0[k]/=sqrt(xxc);
                }
                for (k=0; k<3; k++) {
                    // calculate triangle points at the base of prism
                    pa[k]=vert[i][k]+a0[k]*dedge;
                    pb[k]=vert[i][k]+b0[k]*dedge;
                    pc[k]=vert[i][k]+c0[k]*dedge;
                }
                // set boundaries of suspicious region
                mini[0]=2*xmax; maxi[0]=-2*xmax;
                mini[1]=2*ymax; maxi[1]=-2*ymax;
                mini[2]=2*zmax; maxi[2]=-2*zmax;
                for (k=0; k<3; k++) {
                    if (vert[vinc[i][j]][k]<mini[k]) {
                        mini[k]=vert[vinc[i][j]][k];
                    }
                    if (vert[i][k]<mini[k]) mini[k]=vert[i][k];
                    if (vert[vinc[i][j]][k]>maxi[k]) {
                        maxi[k]=vert[vinc[i][j]][k];
                    }
                    if (vert[i][k]>maxi[k]) maxi[k]=vert[i][k];
                }
                for (k=0; k<3; k++) {
                    mini[k]-=dedge;
                    maxi[k]+=dedge;
                }
                for (k=(int) floor(mini[0]); k<=ceil(maxi[0]); k++)
                    for (l=(int) floor(mini[1]); l<=ceil(maxi[1]); l++)
                        for (m=(int) floor(mini[2]); m<=ceil(maxi[2]); m++){
                            // for all points P in suspicious region
                            if (amat[CONFINEX(k)][CONFINEY(l)]\
                                [CONFINEZ(m)]==0) {
                                da=sqrt(pow(vert[i][0]-k,2)+pow(vert[i][1]-\
                                    l,2)+pow(vert[i][2]-m,2));
                                db=sqrt(pow(k-vert[vinc[i][j]][0],2)+\
                                    pow(l-vert[vinc[i][j]][1],2)+\
                                    pow(m-vert[vinc[i][j]][2],2));
                                if (da<de && db<de) {
                                    // calculate distance of triangle from
                                    // plane of point P (it is the same for
                                    // any point of triangle)
                                    da=abs(n0[0]*pa[0]+n0[1]*pa[1]+n0[2]*\
                                        pa[2]-n0[0]*k-n0[1]*l-n0[2]*m);

                                    pp[0]=k; pp[1]=l; pp[2]=m;
                                    for (n=0; n<3; n++) {
                                        //calculate triangle vertices in
                                        //plane of point P
                                        pap[n]=pa[n]+da*n0[n];
                                        pbp[n]=pb[n]+da*n0[n];
                                        pcp[n]=pc[n]+da*n0[n];
                                    }
                                    if (abs(n0[0]*pap[0]+n0[1]*pap[1]+\
                                        n0[2]*pap[2]-n0[0]*k-n0[1]*l-\
                                        n0[2]*m) > 1e-2) {
                                        for (n=0; n<3; n++) {
                                            //calculate triangle vertices in
                                            //plane of point P
                                            pap[n]=pa[n]-da*n0[n];
                                            pbp[n]=pb[n]-da*n0[n];
                                            pcp[n]=pc[n]-da*n0[n];
                                        }
                                    }
                                    if (abs(n0[0]*pap[0]+n0[1]*pap[1]+\
                                        n0[2]*pap[2]-n0[0]*k-n0[1]*l-\
                                        n0[2]*m) > 1e-2) {
                                        cout << "not in plane" << endl;
                                        exit(1);
                                    }
                                    if (PointInTriangle(pp,pap,pbp,pcp)) {
                                        //make him solid
                                        amat[CONFINEX(k)][CONFINEY(l)]\
                                            [CONFINEZ(m)]=1;
                                    }
                                }
                            }
                        }
            }
        }
    }
    return;
}
