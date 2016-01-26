/*! \file
	\brief Defines global variables, macros, templates and namespace.
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#ifndef GLOBALS_H
#define GLOBALS_H
//! Absolute value
#define ABS(x)		(((x) < 0) ? -(x) : (x))
//! Lesser of two values
#define MIN(x,y)	(((x) < (y)) ? (x) : (y))
//! Periodic boundary conditions in `x`.
#define CONFINEX(x)		(((x) >= 0) ? (x) % nx : (x) + nx)
//! Periodic boundary conditions in `y`.
#define CONFINEY(x)		(((x) >= 0) ? (x) % ny : (x) + ny)
//! Periodic boundary conditions in `z`.
#define CONFINEZ(x)		(((x) >= 0) ? (x) % nz : (x) + nz)
//! Respect periodic boundary conditions in matrix.
#define AMAT(i,j,k)	        amat[CONFINEX(i)][CONFINEY(j)][CONFINEZ(k)]
//! Respect periodic boundary conditions cell volume matrix.
#define V_CELL(m,i,j,k) v_cell[(m)][CONFINEX(i)][CONFINEY(j)][CONFINEZ(k)]
//! signum
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
//! struct for passing paramaters to GSL function
struct fn1_params {int sv; int incmax; int vmax; double **vert; int **vinc; \
    int ***smat; bool report;};
//! namespace with global variables
namespace globals {
    extern bool createNodes; //!< create struts at cell vertices
    extern bool createEdges; //!< create struts at cell edges
    extern bool openCell; //!< make open cell foam
    extern bool save_dat; //!< save file in old dx style
    extern bool save_vtk; //!< save file in new paraview style
    extern double dstrut; //!< parameter influencing size of struts in cell vertices
    extern double dedge; //!< parameter influencing size of struts in cell edges
    extern double strutPorosity; //!< desired porosity of only struts
    //! number between 0.0 and 1.0,
    //! how much are positions of cell centers perturbated from even position
    extern double RANDOM;
    //! grid type for cell centers, 1=cubic grid,
    //! 2=hexagonal lattice, ABAB...,3=hexagonal lattice, ABCABC...
    extern int grid;
    extern int nx; //!< domain size in X
    extern int ny; //!< domain size in Y
    extern int nz; //!< domain size in Z
    extern int sx; //!< size of cells in X
    extern int sy; //!< size of cells in Y
    extern int sz; //!< size of cells in Z
    extern bool save_voro_diag1; //!< gnuplot Voronoi diagram
    extern bool save_voro_diag2; //!< alternative gnuplot Voronoi diagram
    extern bool import_vtk; //!< import morphology from vtk
    extern bool progress_report; //!< show detailed progress report
}
#endif
