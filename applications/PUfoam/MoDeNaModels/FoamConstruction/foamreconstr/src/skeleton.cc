/*! \file
	\brief Functions for geometric tessellation. Vertices and edges are used.
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "voro++.hh"
#include "geometry.hh"
using namespace std;
using namespace voro;
using namespace globals;
//! Make geometric tessellation using `voro++` and store it to a file.
void makeFoamSkeleton(\
    string filename /**< [in] Name of output file with tessellation */,\
    int ncell /**< [in] number of cells */,\
    int *center_x /**< [out] `x` position of seeds */,\
    int *center_y /**< [out] `y` position of seeds */,\
    int *center_z /**< [out] `z` position of seeds */,\
    bool report /**< [in] show output */)
{
    // Uses Voro++ to create foam skeleton - tessellation on provided seeds.
    // Skeleton is saved to gnuplot file.
    int i;
    double xmin=0,xmax=nx;
    double ymin=0,ymax=ny;
    double zmin=0,zmax=nz;
    const int n_x=6,n_y=6,n_z=6; // Set up the number of blocks that
    // the container is divided into (voro++)
    if (report) {
        cout << "voro++ is creating Voronoi tesselation..." << endl;
    }
    // Create a container with the geometry given above, and make it
    // non-periodic in each of the three coordinates. Allocate space for
    // eight particles within each computational block
    container con(\
        xmin,xmax,ymin,ymax,zmin,zmax,n_x,n_y,n_z,true,true,true,8);
    for (i=0; i<ncell; i++) { //put cell centers in container
        con.put(i,center_x[i],center_y[i],center_z[i]);
    }
    // Save the Voronoi network of all the particles to text files
    con.draw_cells_gnuplot(filename.c_str());
}
//! Load geometric tessellation, move it to the domain and store it.
void importFoamSkeleton(\
    string filename /**< [in] Name of output file with tessellation */,\
    double **vert /**< [out] vertex positions */,\
    int **vinc /**< [out] indexes of connected vertices */,\
    int vmax /**< [in] maximum number of vertices */,\
    int incmax /**< [in] maximum number of vertex connections */,\
    int &sv /**< [out] number of vertices */,\
    bool report /**< [in] show output */)
{
    // Imports foam skeleton from gnuplot file. Creates `vert` and `vinc` and
    // calculates `sv`.
    int i,j;
    double xmin=0,xmax=nx;
    double ymin=0,ymax=ny;
    double zmin=0,zmax=nz;
    double x[6],y[6],z[6]; //vertex coordinates
    int vind[4]; //indices of incident vertices
    ifstream fin;
    string line;
    bool found;
    double eps=1e-2; //tolerance for close voronoi vertices
    if (report) {
        cout << "loading cell vertices and edges" << endl;
    }
    fin.open(filename);
    if (!fin.is_open()) {
        cout << "can't open gnuplot file with foam skeleton" << endl;
        exit(1);
    }
    for (i=0; i<vmax; i++) { //initialize vert
        vert[i][0]=0;
        vert[i][1]=0;
        vert[i][2]=0;
        for (j=0; j<incmax; j++) vinc[i][j]=-1;
    }
    sv=0; //number of stored vertices
    getline(fin,line);
    istringstream tmp(line);
    if (!import_vtk) {
        tmp >> x[0]; tmp >> y[0]; tmp >> z[0];
    } else {
        tmp >> z[0]; tmp >> y[0]; tmp >> x[0];
        x[0]=x[0]*nx;
        y[0]=y[0]*ny;
        z[0]=z[0]*nz;
    }
    while (getline(fin,line)) {
    //store coordinates of vertices and their incidence to each other
        if (!line.empty()) {
            istringstream tmp(line);
            if (!import_vtk) {
                tmp >> x[1]; tmp >> y[1]; tmp >> z[1];
            } else {
                tmp >> z[1]; tmp >> y[1]; tmp >> x[1];
                x[1]=x[1]*nx;
                y[1]=y[1]*ny;
                z[1]=z[1]*nz;
            }
            if (in_domain(x[0],y[0],z[0]) && in_domain(x[1],y[1],z[1])) {
            //both vertices are in the domain
                for (i=0; i<2; i++) {
                    found=false;
                    for (j=0; j<sv; j++) {
                        if (abs(x[i]-vert[j][0])<eps && \
                            abs(y[i]-vert[j][1])<eps && \
                            abs(z[i]-vert[j][2])<eps) {
                                found=true;
                                vind[i]=j; //store the vertex index
                                break;
                            }
                    }
                    if (found == false) {
                        //store the vertex coordinates to the end of matrix
                        vert[sv][0]=x[i];
                        vert[sv][1]=y[i];
                        vert[sv][2]=z[i];
                        vind[i]=sv; //store the vertex index
                        sv++;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 1 is incident to vertex 0
                    if (vinc[vind[0]][i] == vind[1]) break;
                    if (vinc[vind[0]][i] == -1) {
                        vinc[vind[0]][i] = vind[1];
                        break;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 0 is incident to vertex 1
                    if (vinc[vind[1]][i] == vind[0]) break;
                    if (vinc[vind[1]][i] == -1) {
                        vinc[vind[1]][i] = vind[0];
                        break;
                    }
                }
            } else if ((in_domain(x[0],y[0],z[0]) && \
                !in_domain(x[1],y[1],z[1])) || (!in_domain(x[0],y[0],z[0]) \
                && in_domain(x[1],y[1],z[1]))) {
                    //only one of the vertices is in the domain
                if (x[0]<0 || x[1]<0) { //mirror the edge over the domain
                    x[2]=x[0]+xmax;
                    x[3]=x[1]+xmax;
                } else if (x[0]>xmax || x[1]>xmax) {
                    x[2]=x[0]-xmax;
                    x[3]=x[1]-xmax;
                } else {
                    x[2]=x[0];
                    x[3]=x[1];
                }
                if (y[0]<0 || y[1]<0) {
                    y[2]=y[0]+ymax;
                    y[3]=y[1]+ymax;
                } else if (y[0]>ymax || y[1]>ymax) {
                    y[2]=y[0]-ymax;
                    y[3]=y[1]-ymax;
                } else {
                    y[2]=y[0];
                    y[3]=y[1];
                }
                if (z[0]<0 || z[1]<0) {
                    z[2]=z[0]+zmax;
                    z[3]=z[1]+zmax;
                } else if (z[0]>zmax || z[1]>zmax) {
                    z[2]=z[0]-zmax;
                    z[3]=z[1]-zmax;
                } else {
                    z[2]=z[0];
                    z[3]=z[1];
                }
                for (i=0; i<4; i++) {
                    found=false;
                    for (j=0; j<sv; j++) {
                        if (abs(x[i]-vert[j][0])<eps && \
                            abs(y[i]-vert[j][1])<eps && \
                            abs(z[i]-vert[j][2])<eps) {
                            found=true;
                            vind[i]=j; //store the vertex index
                            break;
                        }
                    }
                    if (found == false) {
                        //store the vertex coordinates to the end of matrix
                        vert[sv][0]=x[i];
                        vert[sv][1]=y[i];
                        vert[sv][2]=z[i];
                        vind[i]=sv; //store the vertex index
                        sv++;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 1 is incident to vertex 0
                    if (vinc[vind[0]][i] == vind[1]) break;
                    if (vinc[vind[0]][i] == -1) {
                        vinc[vind[0]][i] = vind[1];
                        break;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 0 is incident to vertex 1
                    if (vinc[vind[1]][i] == vind[0]) break;
                    if (vinc[vind[1]][i] == -1) {
                        vinc[vind[1]][i] = vind[0];
                        break;
                    }
                }

                for (i=0; i<incmax; i++) {
                    //vertex 3 is incident to vertex 2
                    if (vinc[vind[2]][i] == vind[3]) break;
                    if (vinc[vind[2]][i] == -1) {
                        vinc[vind[2]][i] = vind[3];
                        break;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 2 is incident to vertex 3
                    if (vinc[vind[3]][i] == vind[2]) break;
                    if (vinc[vind[3]][i] == -1) {
                        vinc[vind[3]][i] = vind[2];
                        break;
                    }
                }
            } else {
                if (x[0]<0) {
                    //mirror the edge over the domain for the vertex 0
                    x[2]=x[0]+xmax;
                    x[3]=x[1]+xmax;
                } else if (x[0]>xmax) {
                    x[2]=x[0]-xmax;
                    x[3]=x[1]-xmax;
                } else {
                    x[2]=x[0];
                    x[3]=x[1];
                }
                if (y[0]<0) {
                    y[2]=y[0]+ymax;
                    y[3]=y[1]+ymax;
                } else if (y[0]>ymax) {
                    y[2]=y[0]-ymax;
                    y[3]=y[1]-ymax;
                } else {
                    y[2]=y[0];
                    y[3]=y[1];
                }
                if (z[0]<0) {
                    z[2]=z[0]+zmax;
                    z[3]=z[1]+zmax;
                } else if (z[0]>zmax) {
                    z[2]=z[0]-zmax;
                    z[3]=z[1]-zmax;
                } else {
                    z[2]=z[0];
                    z[3]=z[1];
                }
                if (x[1]<0) { //mirror the edge over the domain
                    x[4]=x[0]+xmax;
                    x[5]=x[1]+xmax;
                } else if (x[1]>xmax) {
                    x[4]=x[0]-xmax;
                    x[5]=x[1]-xmax;
                } else {
                    x[4]=x[0];
                    x[5]=x[1];
                }
                if (y[1]<0) {
                    y[4]=y[0]+ymax;
                    y[5]=y[1]+ymax;
                } else if (y[1]>ymax) {
                    y[4]=y[0]-ymax;
                    y[5]=y[1]-ymax;
                } else {
                    y[4]=y[0];
                    y[5]=y[1];
                }
                if (z[1]<0) {
                    z[4]=z[0]+zmax;
                    z[5]=z[1]+zmax;
                } else if (z[1]>zmax) {
                    z[4]=z[0]-zmax;
                    z[5]=z[1]-zmax;
                } else {
                    z[4]=z[0];
                    z[5]=z[1];
                }
                for (i=2; i<6; i++) {
                    found=false;
                    for (j=0; j<sv; j++) {
                        if (abs(x[i]-vert[j][0])<eps && \
                            abs(y[i]-vert[j][1])<eps && \
                            abs(z[i]-vert[j][2])<eps) {
                            found=true;
                            vind[i-2]=j; //store the vertex index
                            break;
                        }
                    }
                    if (found == false) {
                        //store the vertex coordinates to the end of matrix
                        vert[sv][0]=x[i];
                        vert[sv][1]=y[i];
                        vert[sv][2]=z[i];
                        vind[i-2]=sv; //store the vertex index
                        sv++;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 1 is incident to vertex 0
                    if (vinc[vind[0]][i] == vind[1]) break;
                    if (vinc[vind[0]][i] == -1) {
                        vinc[vind[0]][i] = vind[1];
                        break;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 0 is incident to vertex 1
                    if (vinc[vind[1]][i] == vind[0]) break;
                    if (vinc[vind[1]][i] == -1) {
                        vinc[vind[1]][i] = vind[0];
                        break;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 3 is incident to vertex 2
                    if (vinc[vind[2]][i] == vind[3]) break;
                    if (vinc[vind[2]][i] == -1) {
                        vinc[vind[2]][i] = vind[3];
                        break;
                    }
                }
                for (i=0; i<incmax; i++) {
                    //vertex 2 is incident to vertex 3
                    if (vinc[vind[3]][i] == vind[2]) break;
                    if (vinc[vind[3]][i] == -1) {
                        vinc[vind[3]][i] = vind[2];
                        break;
                    }
                }
            }
            x[0]=x[1]; y[0]=y[1]; z[0]=z[1];
        } else {
            getline(fin,line); //skip the line,
            //there are always two empty lines together
            getline(fin,line); //refresh first vertex
            istringstream tmp(line);
            if (!import_vtk) {
                tmp >> x[0]; tmp >> y[0]; tmp >> z[0];
            } else {
                tmp >> z[0]; tmp >> y[0]; tmp >> x[0];
                x[0]=x[0]*nx;
                y[0]=y[0]*ny;
                z[0]=z[0]*nz;
            }
        }
    }
    fin.close();
}
