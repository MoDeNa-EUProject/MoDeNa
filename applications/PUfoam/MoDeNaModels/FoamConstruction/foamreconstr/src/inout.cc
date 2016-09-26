/*! \file
	\brief Functions for input and output operations.
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#include "globals.hh"
#include <fstream>
#include <iostream>
#include <string.h>
#include "allocation.hh"
#include "geometry.hh"
using namespace std;
using namespace globals;
//! Reads file with options and parameters
//!
//! Reads parameters from file and stores them to global variables. These
//! parameters determine what will be done.
void readParameters(\
    //! [in] Name of input file with parameters that determine what will be done
    string filename,\
    //! [in,out] Name of output file with morphology without extension
    string &outputFilename,\
    //! [in,out] Name of input VTK file with morphology
    string &VTKInputFilename,\
    //! [in,out] Name of input/output file with tessellation - gnuplot diagram
    string &GnuplotSkeletonFilename,\
    //! [in,out] Name of output file with tessellation - alt. gnuplot diagram
    string &GnuplotAltSkeletonFilename,\
    //! [in,out] Name of output file with morphology descriptors
    string &descriptorsFilename,\
    //! [in,out] Name of output file with parameters
    string &parametersFilename)
{
    // Reads inputs from text file. Returns filenames. Rest is put into global
    // variables.
    ifstream fin;
    fin.open(filename);
        if (!fin.is_open()) {
            cout << "Can't open file " << filename << endl;
            exit(1);
        }
        fin >> createNodes; fin.ignore(256,'\n');
        fin >> createEdges; fin.ignore(256,'\n');
        fin >> openCell; fin.ignore(256,'\n');
        fin >> save_dat; fin.ignore(256,'\n');
        fin >> save_vtk; fin.ignore(256,'\n');
        fin >> dstrut; fin.ignore(256,'\n');
        fin >> dedge; fin.ignore(256,'\n');
        fin >> strutPorosity; fin.ignore(256,'\n');
        fin >> RANDOM; fin.ignore(256,'\n');
        fin >> grid; fin.ignore(256,'\n');
        fin >> nx; fin.ignore(256,'\n');
        fin >> ny; fin.ignore(256,'\n');
        fin >> nz; fin.ignore(256,'\n');
        fin >> sx; fin.ignore(256,'\n');
        fin >> sy; fin.ignore(256,'\n');
        fin >> sz; fin.ignore(256,'\n');
        fin >> save_voro_diag1; fin.ignore(256,'\n');
        fin >> save_voro_diag2; fin.ignore(256,'\n');
		fin >> import_vtk; fin.ignore(256,'\n');
        fin >> progress_report; fin.ignore(256,'\n');
        fin >> outputFilename; fin.ignore(256,'\n');
        fin >> VTKInputFilename; fin.ignore(256,'\n');
        fin >> GnuplotSkeletonFilename; fin.ignore(256,'\n');
        fin >> GnuplotAltSkeletonFilename; fin.ignore(256,'\n');
        fin >> descriptorsFilename; fin.ignore(256,'\n');
        fin >> parametersFilename; fin.ignore(256,'\n');
    fin.close();
}
//! Allocates the matrix from dimensions read from VTK file
int ***allocateFromVTK(\
    string filename /**< [in] Name of input VTK file with morphology */,\
    int ***amat /**< [in,out] matrix */,\
    bool report /**< [in] show output */)
{
    int i,j,k;
    ifstream fin;
    string line;
    // read vtk file
    if (report) {
        cout << "reading vtk file" << endl;
    }
    fin.open(filename);
    // read header
    int vtkx,vtky,vtkz;
    for (i=0;i<9;i++) {
        getline(fin, line);
        if (i==4) {
            char word[20];
            sscanf(line.c_str(),"%s %d %d %d",word,&vtkx,&vtky,&vtkz);
        }
    }
    if (report) {
        cout << "dimensions: " << vtkx << ", " << vtky << ", " << vtkz << endl;
        cout << "allocating" << endl;
    }
    amat = alloc_3Dmatrix (vtkx, vtky, vtkz);
    for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++)
        amat[i][j][k] = 0;
    nx=vtkx;
    ny=vtky;
    nz=vtkz;
    return amat;
}
//! Reads morphology from VTK file. Changes `amat`.
void importFromVTK(\
    string filename /**< [in] Name of input VTK file with morphology */,\
    int ***amat /**< [in,out] matrix */,\
    bool report /**< [in] show output */)
{
    int i,j,k,l;
    ifstream fin;
    string line;
    int *bmat;
    float a,b,c;
    // read vtk file
    if (report) {
        cout << "reading vtk file" << endl;
    }
    fin.open(filename);
    // read header
    int vtkx,vtky,vtkz;
    for (i=0;i<9;i++) {
        getline(fin, line);
        if (i==4) {
            char word[20];
            sscanf(line.c_str(),"%s %d %d %d",word,&vtkx,&vtky,&vtkz);
        }
    }
    if (report) {
        cout << "dimensions: " << vtkx << ", " << vtky << ", " << vtkz << endl;
        cout << "loading data file" << endl;
    }
    bmat = (int *)calloc((size_t)vtkx*vtky*vtkz, sizeof(int));
    // read data
    getline(fin, line);
    sscanf(line.c_str(),"%e %e %e",&a,&b,&c);
    // 3 values on first data line
    bmat[0]=int(a+0.5);
    bmat[1]=int(b+0.5);
    bmat[2]=int(c+0.5);
    l=3;
    for (i=1;i<vtkx*vtky*vtkz/2-1;i++) {
        // 2 values on next lines
        getline(fin, line);
        sscanf(line.c_str(),"%e %e",&a,&b);
        bmat[2*i+1]=int(a+0.5);
        bmat[2*i+2]=int(b+0.5);
        l=l+2;
    }
    //check that we read the right amount of values
    // cout << vtkx*vtky*vtkz << " " << l << endl;
    // 1 value on last line
    getline(fin, line);
    sscanf(line.c_str(),"%e",&a);
    bmat[vtkx*vtky*vtkz-1]=int(a+0.5);
    fin.close();
    // make 3D array from 1D array
    l=0;
    for (i=0;i<vtkz;i++) {
        for (j=0;j<vtky;j++) {
            for (k=0;k<vtkx;k++) {
                if (amat[i][j][k]==0) {
                    amat[i][j][k]=bmat[l];
                }
                l++;
            }
        }
    }
    free(bmat);
    nx=vtkx;
    ny=vtky;
    nz=vtkz;
}
//! Saves morphology to VTK file.
void saveToVTK(\
    const char* filename /**< [in] Name of output file with morphology */,\
    int ***amat /**< [in] matrix */,\
    bool report /**< [in] show output */)
{
    int i,j,k;
    FILE *strmo;
    // Output to file in Paraview style.
    if (report) {
        cout << "saving in Paraview style..." << endl;
    }
    strmo = fopen(filename, "w");
    if (strmo == NULL) {
        fprintf (stderr, "Can't open file %s\n",filename);
        exit(13);
    }

    fprintf (strmo, "# vtk DataFile Version 3.0\n");
    fprintf (strmo, "vtkfile\n");
    fprintf (strmo, "ASCII\n");
    fprintf (strmo, "DATASET STRUCTURED_POINTS\n");
    fprintf (strmo, "DIMENSIONS %d %d %d\n",nx,ny,nz);
    fprintf (strmo, "ORIGIN 1 1 1\n");
    fprintf (strmo, "SPACING %g %g %g\n",1.0/nx,1.0/ny,1.0/nz);
    fprintf (strmo, "POINT_DATA %d\n",nx*ny*nz);
    fprintf (strmo, "SCALARS values int\n");
    fprintf (strmo, "LOOKUP_TABLE default\n");

    for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++)
        fprintf (strmo, " %c", (amat[i][j][k] >= 1) ? '0' : '1');
        //fprintf (strmo, " %d", amat[i][j][k]);
    fprintf (strmo, "\n");
    }
    fprintf (strmo, "\n");
    }
    fprintf (strmo, "\n");

    fclose (strmo);
}
//! Saves morphology to DX file.
void saveToDX(\
    const char* filename /**< [in] Name of output file with morphology */,\
    int ***amat /**< [in] matrix */,\
    bool report /**< [in] show output */)
{
    int i,j,k;
    FILE *strmo;
    // Output to file in DX style.
    if (report) {
        cout << "saving in DX style..." << endl;
    }
    strmo = fopen(filename, "w");
    if (strmo == NULL) {
        fprintf (stderr, "Can't open file %s\n",filename);
        exit(13);
    }

    fprintf (strmo, "NCX =    %d    ;\n", nx);
    fprintf (strmo, "NCY =    %d    ;\n", ny);
    fprintf (strmo, "NCZ =    %d    ;\n", nz);

    fprintf (strmo, "A = [\n");
    for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
    for (k = 0; k < nz; k++)
        fprintf (strmo, " %c", (amat[i][j][k] >= 1) ? '0' : '1');
        //fprintf (strmo, " %d", amat[i][j][k]);
    fprintf (strmo, "\n");
    }
    fprintf (strmo, "\n");
    }
    fprintf (strmo, "];\n");

    fclose (strmo);
}
//! Saves tessellation diagram to gnuplot file.
void saveToGnuplot(\
    string filename /**< [in] Name of output file with tessellation */,\
    int sv /**< [in] number of vertices */,\
    int incmax /**< [in] maximum number of connections to one vertex */,\
    double **vert /**< [in] vertex positions */,\
    int **vinc /**< [in] indexes of connected vertices */,\
    bool report /**< [in] show output */)
{
    // Stores cell edges incident to voronoi vertices in the domain
    ofstream fout;
    int i,j;
    if (report) {
        cout << "saving shifted Voronoi diagram" << endl;
    }
    fout.open(filename);
    if (!fout.is_open()) {
        cout << "can't open file " << filename << endl;
        exit(1);
    }
    for (i=0; i<sv; i++) {
        if (in_domain(vert[i][0],vert[i][1],vert[i][2])) {
            for (j=0; j<incmax; j++) {
                if (vinc[i][j] != -1) {
                    fout << vert[i][0] << " " << vert[i][1] << " " << \
                        vert[i][2] << endl;
                    fout << vert[vinc[i][j]][0] << " " << \
                        vert[vinc[i][j]][1] << " " << \
                        vert[vinc[i][j]][2] << endl;
                    fout << endl;
                    fout << endl;
                } else {
                    break;
                }
            }
        }
    }
    fout.close();
}
//! Saves morphology descriptors (porosity and strut content) to a file.
void saveDescriptors(\
    string filename /**< [in] Output file with morphology descriptors */,\
    double por /**< [in] porosity */,\
    double fs /**< [in] strut content */,\
    bool report /**< [in] show output */)
{
    ofstream fout;
    if (report) {
        cout << "saving descriptors" << endl;
    }
    fout.open(filename);
    if (!fout.is_open()) {
        cout << "can't open file " << filename << endl;
        exit(1);
    }
    fout << por << endl;
    fout << fs << endl;
    fout.close();
}
//! Saves parameters (`dedge`) to a file.
void saveParameters(\
    string filename /**< [in] Name of output file with parameters */,\
    double dedge /**< [in] edge strut size parameter */,\
    bool report /**< [in] show output */)
{
    ofstream fout;
    if (report) {
        cout << "saving parameters" << endl;
    }
    fout.open(filename);
    if (!fout.is_open()) {
        cout << "can't open file " << filename << endl;
        exit(1);
    }
    fout << dedge << endl;
    fout.close();
}
