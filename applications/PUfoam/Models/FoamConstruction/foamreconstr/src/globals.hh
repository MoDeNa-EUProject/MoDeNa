#ifndef GLOBALS_H
#define GLOBALS_H
//
// Macros.
//
#define ABS(x)		(((x) < 0) ? -(x) : (x))
#define MIN(x,y)	(((x) < (y)) ? (x) : (y))
//
// Periodic boundary conditions are done by these simple macros.
//
#define CONFINEX(x)		(((x) >= 0) ? (x) % nx : (x) + nx)
#define CONFINEY(x)		(((x) >= 0) ? (x) % ny : (x) + ny)
#define CONFINEZ(x)		(((x) >= 0) ? (x) % nz : (x) + nz)
#define AMAT(i,j,k)	        amat[CONFINEX(i)][CONFINEY(j)][CONFINEZ(k)]
#define V_CELL(m,i,j,k) v_cell[(m)][CONFINEX(i)][CONFINEY(j)][CONFINEZ(k)]
//
// Templates.
//
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
namespace globals {
    extern bool createNodes; //create struts at cell vertices
    extern bool createEdges; //create struts at cell edges
    extern bool openCell; //make open cell foam
    extern bool save_dat; //save file in old dx style
    extern bool save_vtk; //save file in new paraview style
    extern float dstrut; //parameter influencing size of struts in cell vertices
    extern float dedge; //parameter influencing size of struts in cell edges
    extern float RANDOM; //number between 0.0 and 1.0
    //how much are positions of cell centers perturbated from even position
    extern int grid; //grid type for cell centers, 1=cubic grid,
    //2=hexagonal lattice, ABAB...,3=hexagonal lattice, ABCABC...
    extern int nx; //domain size in X
    extern int ny; //domain size in Y
    extern int nz; //domain size in Z
    extern int sx; //size of cells in X
    extern int sy; //size of cells in Y
    extern int sz; //size of cells in Z
    extern bool save_voro_diag1; //gnuplot Voronoi diagram
    extern bool save_voro_diag2; //alternative gnuplot Voronoi diagram
    extern bool import_vtk; //import morphology from vtk
}
#endif
