#include "globals.hh"
#include <fstream>
#include <iostream>
#include <string.h>
#include "allocation.hh"
#include "geometry.hh"
using namespace std;
using namespace globals;
void readParameters(string filename, string &outputFilename, \
    string &VTKInputFilename, string &GnuplotSkeletonFilename, \
    string &GnuplotAltSkeletonFilename, string &descriptorsFilename) {
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
        fin >> outputFilename; fin.ignore(256,'\n');
        fin >> VTKInputFilename; fin.ignore(256,'\n');
        fin >> GnuplotSkeletonFilename; fin.ignore(256,'\n');
        fin >> GnuplotAltSkeletonFilename; fin.ignore(256,'\n');
        fin >> descriptorsFilename; fin.ignore(256,'\n');
    fin.close();
}
int ***allocateFromVTK(string filename, int ***amat) {
    int i,j,k;
    ifstream fin;
    string line;
    // read vtk file
    cout << "reading vtk file" << endl;
    fin.open(filename);
    // read header
    int vtkx,vtky,vtkz;
    for (i=0;i<9;i++) {
        getline(fin, line);
        cout << line << endl;
        if (i==4) {
            char word[20];
            sscanf(line.c_str(),"%s %d %d %d",word,&vtkx,&vtky,&vtkz);
        }
    }
    cout << "dimensions: " << vtkx << ", " << vtky << ", " << vtkz << endl;
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
void importFromVTK(string filename, int ***amat) {
    int i,j,k,l;
    ifstream fin;
    string line;
    int *bmat;
    float a,b,c;
    // read vtk file
    cout << "reading vtk file" << endl;
    fin.open(filename);
    // read header
    int vtkx,vtky,vtkz;
    for (i=0;i<9;i++) {
        getline(fin, line);
        cout << line << endl;
        if (i==4) {
            char word[20];
            sscanf(line.c_str(),"%s %d %d %d",word,&vtkx,&vtky,&vtkz);
        }
    }
    cout << "dimensions: " << vtkx << ", " << vtky << ", " << vtkz << endl;
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
    cout << vtkx*vtky*vtkz << " " << l << endl;
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
void saveToVTK(const char* filename, int ***amat) {
    int i,j,k;
    FILE *strmo;
    // Output to file in Paraview style.
    cout << "saving in Paraview style..." << endl;
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
    fprintf (strmo, "ORIGIN 0 0 0\n");
    fprintf (strmo, "SPACING 1 1 1\n");
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
void saveToDX(const char* filename, int ***amat) {
    int i,j,k;
    FILE *strmo;
    // Output to file in DX style.
    cout << "saving in DX style..." << endl;
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
void saveToGnuplot(string filename, int sv, int incmax, float **vert, \
                   int **vinc) {
    // Stores cell edges incident to voronoi vertices in the domain
    ofstream fout;
    int i,j;
    cout << "saving shifted Voronoi diagram..." << endl;
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
void saveDescriptors(string filename, double por, double fs) {
    ofstream fout;
    fout.open(filename);
    if (!fout.is_open()) {
        cout << "can't open file " << filename << endl;
        exit(1);
    }
    fout << por << endl;
    fout << fs << endl;
    fout.close();
}
