/*! \file
	\brief Main reconstruction algorithm

	Version: v0.3
	\author Pavel Ferkl
	\author Juraj Kosek
	\ingroup foam_constr
*/
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
//! Reads parameters. Creates struts and walls. Saves foam morphology to a file.
int main(int argc,char *argv[])
{
	int ncell; // number of seeds for tessellation
	int i, j, k;
	int ***amat,***smat; // 3D matrix of voxels
    int *center_x, *center_y, *center_z; // position of seeds
	int vmax=10000; // maximum number of vertices
	int incmax=8; // maximum number of incident vertices to a vertex
	double **vert; // vertex coordinates
	int **vinc; // incident vertices for each vertex
	int sv=0; // final number of vertices
	double por=0; // porosity
	double por_s=0; // porosity of struts only
	double fs=0; // strut content
	int gsl_status; // status return by GSL function
	double mini,maxi; // bounds of optimization interval
	int maxit=3; // maximum number of restarts of optimization
	string inputsFilename="foamreconstr.in"; // Name of input file with
		// parameters that determine what will be done
	string outputFilename; // Name of output file with morphology
		// without extension
	string VTKInputFilename; // Name of input VTK file with morphology
	string GnuplotSkeletonFilename; // Name of input/output file with
		// tessellation - gnuplot diagram
	string GnuplotAltSkeletonFilename; // Name of output file with
		// tessellation - alternative gnuplot diagram
	string descriptorsFilename; // Name of output file with
		// morphology descriptors
	string parametersFilename; // Name of output file with parameters

    // read inputs and save them to global variables
	readParameters(inputsFilename,outputFilename,VTKInputFilename,\
		GnuplotSkeletonFilename,GnuplotAltSkeletonFilename,descriptorsFilename,\
		parametersFilename);
	// define names of output files with morphology
	string VTKOutputFilename=outputFilename+".vtk";
	string DXOutputFilename=outputFilename+".dat";
	if (import_vtk) {
		// only allocate the 3D matrix
		amat=allocateFromVTK(VTKInputFilename,amat,progress_report);
	} else {
		// allocate the 3D matrix and make seeds for the tessellation
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
		createSeeds(ncell,center_x,center_y,center_z,progress_report);
	}
    if (createNodes || createEdges) {
		if (!import_vtk) {
		    // geometrical voronoi tesselation by voro++
			// creates gnuplot diagram with the tessellation
			makeFoamSkeleton(GnuplotSkeletonFilename,ncell,center_x,center_y,\
				center_z,progress_report);
		}
		vert=alloc_dmatrix(vmax,3);
		vinc=alloc_matrix(vmax,incmax);
		// read the gnuplot diagram with the tessellation
		// make the alternative diagram with the tessellation
		importFoamSkeleton(GnuplotSkeletonFilename,vert,vinc,vmax,incmax,sv,\
			progress_report);
        if (!save_voro_diag1) {
			// delete the gnuplot diagram with the tessellation from disk
            remove(GnuplotSkeletonFilename.c_str());
        }
        // save alternative gnuplot diagram of the tesselation
        if (save_voro_diag2) {
			saveToGnuplot(GnuplotAltSkeletonFilename,sv,incmax,vert,vinc,\
				progress_report);
        }
		// allocate matrix for struts
	    smat = alloc_3Dmatrix (nx, ny, nz);
	    if (smat == NULL) {
			fprintf (stderr, "Insufficient memory.\n");
			return 9;
	    }
		// struct for passing variables to GSL function
		struct fn1_params params = {sv,incmax,vmax,vert,vinc,smat,\
			progress_report};
		// initial interval where we are looking for optimal dedge
		if (dedge == 0) {
			dedge=2;
		}
		mini=0.5*dedge;
		maxi=1.5*dedge;
		for (i=0;i<maxit;i++) {
			// finds dedge, which gives desired porosity of struts
			gsl_status = optim(&params,dedge,mini,maxi,progress_report);
			if (gsl_status==0) {
				break;
			} else {
				mini=0.5*mini;
				maxi=1.5*maxi;
				if (i==maxit-1) {
					cout << "Didn't find initial interval" << endl;
					exit(1);
				}
			}
		}
		// save dedge so that we can use it as initial guess next time
		saveParameters(parametersFilename,dedge,progress_report);
    }
    por_s=porosity(smat);
	if (progress_report) {
    	cout << "porosity of struts only " << por_s << endl;
	}
	vert=free_dmatrix(vert);
	vinc=free_matrix(vinc);
	if (!openCell) {
		// create walls
		if (import_vtk) {
			// walls from file
			importFromVTK(VTKInputFilename,amat,progress_report);
		} else {
			// walls from Voronoi tessellation
			makeWalls(amat,ncell,center_x,center_y,center_z,progress_report);
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
	fs=(1-por_s)/(1-por);
	if (progress_report) {
	    cout << "porosity " << por << endl;
	    cout << "polymer in struts " << fs << endl;
	}
	// save and exit
	saveDescriptors(descriptorsFilename,por,fs,progress_report);
    if (save_dat) {
        saveToDX(DXOutputFilename.c_str(),smat,progress_report);
    }
    if (save_vtk) {
		saveToVTK(VTKOutputFilename.c_str(),smat,progress_report);
    }
	amat=free_3Dmatrix(amat);
	amat=free_3Dmatrix(smat);
	exit(0);
}
