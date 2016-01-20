/*! \file
	\brief Functions for the initialization of seeds.
	\author Pavel Ferkl
    \author Juraj Kosek
	\ingroup foam_constr
*/
#include "globals.hh"
#include <stdlib.h>
#include <stdio.h>
using namespace globals;
//! Initialization of seeds. Stores positions of seeds.
void createSeeds(\
    int &ncell /**< [in,out] in: max number of cells, out: actual number */,\
    int *center_x /**< [out] `x` position of seeds */,\
    int *center_y /**< [out] `y` position of seeds */,\
    int *center_z /**< [out] `z` position of seeds */,\
    bool report /**< [in] show output */)
{
    // Creates seeds for Voronoi tessellation at centers of cells at cubic or
    // hexagonal grids. Seed positions are store in `center_x`, `center_y` and
    // `center_z`. Updates `ncell`.
    const int GRID_CUBIC=1,GRID_HEXAG_AB=2,GRID_HEXAG_ABC=3;
    int i,j,k,m;
    int ttt;
    if (report) {
        if        (grid == GRID_CUBIC     ) {
            printf("Cubic grid.\n");
        } else if (grid == GRID_HEXAG_AB  ) {
            printf("Hexagonal grid ABABAB...\n");
        } else if (grid == GRID_HEXAG_ABC ) {
            printf("Hexagonal grid ABCABCABC...\n");
        }
        printf("Matrix of %d x %d x %d voxels.\n", nx, ny, nz);
        printf("Spacing of the grid is %d, %d, %d\n", sx, sy, sz);
        printf("Maximum number of cells, ncell = %d\n", ncell);
    }
    if        (grid == GRID_CUBIC     ) {
        m = 0;
        for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
            if ((i % sx == 0) && (j % sy == 0) && \
                (k % sz == 0)) {
                    center_x[m] = i;
                    center_y[m] = j;
                    center_z[m] = k;
                    m++;
                }
    } else if (grid == GRID_HEXAG_AB  ) {
        m = 0;
        for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
            if        ((i %    sx  == 0) &&
                       (j % (2*sy) == 0) &&
                       (k % (2*sz) == 0)	 ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == (sx/2)) &&
                       (j % (2*sy) ==  sy   ) &&
                       (k % (2*sz) ==  0          )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == (sx/2)  ) &&
                       (j % (2*sy) == 1 + sy/3) &&
                       (k % (2*sz) == sz      )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == 0                         ) &&
                       (j % (2*sy) == sy + 1  + sy/3) &&
                       (k % (2*sz) == sz                  )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            }
    } else if (grid == GRID_HEXAG_ABC ) {
        m = 0;
        for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        for (k = 0; k < nz; k++)
            if        ((i %    sx  == 0) &&
                       (j % (2*sy) == 0) &&
                       (k % (3*sz) == 0)	 ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == (sx/2)) &&
                       (j % (2*sy) ==  sy   ) &&
                       (k % (3*sz) ==  0          )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == (sx/2)  ) &&
                       (j % (2*sy) == 1 + sy/3) &&
                       (k % (3*sz) == sz      )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == 0                         ) &&
                       (j % (2*sy) == sy + 1  + sy/3) &&
                       (k % (3*sz) == sz                  )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == 0                 ) &&
                       (j % (2*sy) == 1 + (2*sy)/3) &&
                       (k % (3*sz) == (2*sz)      )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            } else if ((i %    sx  == (sx/2)                 ) &&
                       (j % (2*sy) == sy + 1 + (2*sy)/3) &&
                       (k % (3*sz) == (2*sz)                 )   ) {
                center_x[m] = i;
                center_y[m] = j;
                center_z[m] = k;
                m++;
            }
    }
    if (report) {
        printf("Number of initialized centers of cells = %d\n", m);
    }
    ncell = m;

    /*
     * Random perturbations of the grid.
     */
    for (m = 0; m < ncell; m++) {

        /*
         * center_x[m]
         */
        ttt = rand();
        if (ttt < RAND_MAX/2) {
            center_x[m] += (int)(RANDOM *
                                 (double)ttt / (double)RAND_MAX *
                                 (double)sx);
        } else {
            ttt = ttt - RAND_MAX/2;
            center_x[m] -= (int)(RANDOM *
                                 (double)ttt / (double)RAND_MAX *
                                 (double)sx);
        }
        center_x[m] = CONFINEX( center_x[m] );

        /*
         * center_y[m]
         */
        ttt = rand();
        if (ttt < RAND_MAX/2) {
            center_y[m] += (int)(RANDOM *
                                 (double)ttt / (double)RAND_MAX *
                                 (double)sy);
        } else {
            ttt = ttt - RAND_MAX/2;
            center_y[m] -= (int)(RANDOM *
                                 (double)ttt / (double)RAND_MAX *
                                 (double)sy);
        }
        center_y[m] = CONFINEY( center_y[m] );

        /*
         * center_z[m]
         */
        ttt = rand();
        if (ttt < RAND_MAX/2) {
            center_z[m] += (int)(RANDOM *
                                 (double)ttt / (double)RAND_MAX *
                                 (double)sz);
        } else {
            ttt = ttt - RAND_MAX/2;
            center_z[m] -= (int)(RANDOM *
                                 (double)ttt / (double)RAND_MAX *
                                 (double)sz);
        }
        center_z[m] = CONFINEZ( center_z[m] );

    }
}
