// Foam reconstruction algorithm
// author: Pavel Ferkl
// author: Juraj Kosek
// version: v0.3
// TODO documentation
#include "globals.hh"
#include <iostream>
#include <string.h>
#include "inout.hh"
#include "allocation.hh"
#include "geometry.hh"
#include "skeleton.hh"
#include "walls.hh"
#include "seeds.hh"
#include "struts.hh"
using namespace std;
using namespace globals;

int main(int argc,char *argv[])
{
	int ncell; //number of seeds for tessellation
	int i, j, k;
	int ***amat,***smat; //3D matrix of voxels
    int *center_x, *center_y, *center_z; //position of seeds
	int vmax=10000; //maximum number of vertices
	int incmax=8; //maximum number of incident vertices to a vertex
	float **vert; //vertex coordinates
	int **vinc; //incident vertices for each vertex
	int sv=0; //final number of vertices
	double por=0; //porosity
	double por_s=0; //porosity of struts only
	double fs=0; //strut content
	int gsl_status;
	string inputsFilename="foamreconstr.in";
	string outputFilename;
	string VTKInputFilename;
	string GnuplotSkeletonFilename;
	string GnuplotAltSkeletonFilename;
	string descriptorsFilename;

    // read inputs and save them to global variables
	readParameters(inputsFilename,outputFilename,VTKInputFilename,\
		GnuplotSkeletonFilename,GnuplotAltSkeletonFilename,descriptorsFilename);
	string VTKOutputFilename=outputFilename+".vtk";
	string DXOutputFilename=outputFilename+".dat";
	if (import_vtk) {
		amat=allocateFromVTK(VTKInputFilename,amat);
	} else {
	    ncell =          1 + (nx - 1) / sx ;
	    ncell = ncell * (1 + (ny - 1) / sy);
	    ncell = ncell * (1 + (nz - 1) / sz);
	    /*
	     * Allocate memory:
	     * center_x, center_y, center_z = distance from centers of cells.
	     */
	    center_x =    (int *)calloc((size_t)ncell, sizeof(int));
	    center_y =    (int *)calloc((size_t)ncell, sizeof(int));
	    center_z =    (int *)calloc((size_t)ncell, sizeof(int));
	    if (center_x == NULL || center_y == NULL || center_z == NULL) {
	        fprintf (stderr, "Insufficient memory.\n");
	        return 9;
	    }
	    /*
	     * Allocate the 3D matrix of voxels.
	     * voxel = 0 (pore phase).
	     * voxel > 1 (solid phase).
	     */
	    amat = alloc_3Dmatrix (nx, ny, nz);
	    if (amat == NULL) {
	        fprintf (stderr, "Insufficient memory.\n");
	        return 9;
	    }
	    for (i = 0; i < nx; i++)
	    for (j = 0; j < ny; j++)
	    for (k = 0; k < nz; k++)
	        amat[i][j][k] = 0;
		createSeeds(ncell,center_x,center_y,center_z);
	}
    if (createNodes || createEdges) {
		if (!import_vtk) {
		    // geometrical voronoi tesselation by voro++
			makeFoamSkeleton(GnuplotSkeletonFilename,ncell,center_x,center_y,\
				center_z);
		}
		vert=alloc_fmatrix(vmax,3);
		vinc=alloc_matrix(vmax,incmax);
		importFoamSkeleton(GnuplotSkeletonFilename,vert,vinc,vmax,incmax,sv);
        if (!save_voro_diag1) {
            remove(GnuplotSkeletonFilename.c_str());
        }
        // save alternative gnuplot image of voronoi tesselation
        if (save_voro_diag2) {
			saveToGnuplot(GnuplotAltSkeletonFilename,sv,incmax,vert,vinc);
        }
		// allocate matrix for struts
	    smat = alloc_3Dmatrix (nx, ny, nz);
	    if (smat == NULL) {
			fprintf (stderr, "Insufficient memory.\n");
			return 9;
	    }
		struct fn1_params params = {sv,incmax,vmax,vert,vinc,smat};
		// por=fn1(por,&params);
		gsl_status = optim(&params,dedge);
    }
    por_s=porosity(smat);
    cout << "porosity of struts only " << por_s << endl;
	vert=free_fmatrix(vert);
	vinc=free_matrix(vinc);
	if (!openCell) {
		if (import_vtk) {
			importFromVTK(VTKInputFilename,amat);
		} else {
			makeWalls(amat,ncell,center_x,center_y,center_z);
			free(center_x);
			free(center_y);
			free(center_z);
		}
	}
	for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	for (k = 0; k < nz; k++)
		smat[i][j][k] = smat[i][j][k] + amat[i][j][k];
	por=porosity(smat);
    cout << "porosity " << por << endl;
	fs=(1-por_s)/(1-por);
    cout << "polymer in struts " << fs << endl;
	saveDescriptors(descriptorsFilename,por,fs);
    if (save_dat) {
        saveToDX(DXOutputFilename.c_str(),smat);
    }
    if (save_vtk) {
		saveToVTK(VTKOutputFilename.c_str(),smat);
    }
	amat=free_3Dmatrix(amat);
	amat=free_3Dmatrix(smat);
	exit(0);
}
