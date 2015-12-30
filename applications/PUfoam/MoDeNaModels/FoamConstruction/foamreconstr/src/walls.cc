/*! \file
	\brief Functions for creation of walls from Voronoi tessellation.
	\author Pavel Ferkl
	\author Juraj Kosek
	\ingroup foam_constr
*/
#include "globals.hh"
#include <stdio.h>
#include <stdlib.h>
#include "allocation.hh"
using namespace globals;
//! Some voxels have the same distance from several 'M'
//! centers of cells. Calculate the number 'M'.
#define COUNT(i,j,k,a)	(((V_CELL(0,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(1,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(2,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(3,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(4,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(5,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(6,(i),(j),(k)) == (a)) ? 1 : 0) + \
						 ((V_CELL(7,(i),(j),(k)) == (a)) ? 1 : 0)    )
//! Creates walls based on position of seeds and Voronoi tessellation.
void makeWalls(\
	int ***amat /**< [in,out] voxel matrix */,\
	int ncell /**< [in] number of cells */,\
	int *center_x /**< [out] `x` position of seeds */,\
    int *center_y /**< [out] `y` position of seeds */,\
    int *center_z /**< [out] `z` position of seeds */,\
    bool report /**< [in] show output */)
{
    // Make cell walls based on position of Voronoi seeds. Update `amat`.
    int i,j,k,m;
    int min_distance, m_cell;
    int dx, dy, dz, mx, my, mz;
	int ****v_cell, *v;
	if (report) {
		printf("creating walls\n");
	}
	// v = working array, distances from centers of cells.
	v        =    (int *)calloc((size_t)ncell, sizeof(int));
	/*
	 * Allocate Voronoi set v_cell[8][i][j][k].
	 */
	v_cell = (int ****)calloc((size_t)8, sizeof(int ***));
	for (m = 0; m < 8; m++) {
		v_cell[m] = alloc_3Dmatrix (nx, ny, nz);
		if (v_cell[m] == NULL) {
			fprintf (stderr, "Insufficient memory.\n");
			exit(9);
		}
	}
    /*
     * Initialize the voronoi set v_cell[8][i][j][k].
     */
    for (m = 0; m < 8 ; m++)
    for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++)
        v_cell[m][i][j][k] = -1;

    /*
     * Set the Voronoi set v_cell[8][i][j][k].
     */
    for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++) {
        /*
         * Store the distance of [i,j,k] from all centers of cells
         * in vector v[m].
         */
        for (m = 0; m < ncell; m++) {
            mx = center_x[m];
            my = center_y[m];
            mz = center_z[m];
            dx = MIN(ABS(i - mx), nx - ABS(i - mx));
            dy = MIN(ABS(j - my), ny - ABS(j - my));
            dz = MIN(ABS(k - mz), nz - ABS(k - mz));
            v[m] = dx * dx + dy * dy + dz * dz;
        }
        /*
         * Find the shortest distance from all centers of cells.
         */
        min_distance = 3 * nx * ny * nz;
        for (m = 0; m < ncell; m++)
            if (v[m] < min_distance)
                min_distance = v[m];
        /*
         * Set the Voronoi set v_cell[8][i][j][k] .
         */
        m_cell = 0;
        for (m = 0; m < ncell; m++)
            if (v[m] == min_distance) {
                v_cell[m_cell][i][j][k] = m;
                m_cell++;
                if (m_cell == 8)
                    break;
            }
    }

    /*
     * Set all elements with equal distance from two or more
     * centers of cells to be solid phase.
     */
    m_cell = 0;
    for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++)
        if (v_cell[1][i][j][k] != -1) {
            amat[i][j][k] = 8 - COUNT(i,j,k,-1);
            m_cell++;
        }
	if (report) {
    	printf("Voxels with multiplicity: m_cell = %d\n", m_cell);
	}

    /*
     * Set all elements where at least one of the 6 nearest
     * neighbors is from another Voronoi set to be also wall
     * elements.
     */

    for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
    for (k = 0; k < nz; k++) {
        if		  (amat[i][j][k] > 1) 				       {
            continue;
        } else if (AMAT(i+1,j,k) > 1  				   ||
                   V_CELL(0,i,j,k) != V_CELL(0,i+1,j,k)  ) {
            amat[i][j][k] = 1;
        } else if (AMAT(i-1,j,k) > 1				   ||
                   V_CELL(0,i,j,k) != V_CELL(0,i-1,j,k)  ) {
            amat[i][j][k] = 1;
        } else if (AMAT(i,j+1,k) > 1				   ||
                   V_CELL(0,i,j,k) != V_CELL(0,i,j+1,k)  ) {
            amat[i][j][k] = 1;
        } else if (AMAT(i,j-1,k) > 1				   ||
                   V_CELL(0,i,j,k) != V_CELL(0,i,j-1,k)  ) {
            amat[i][j][k] = 1;
        } else if (AMAT(i,j,k+1) > 1				   ||
                   V_CELL(0,i,j,k) != V_CELL(0,i,j,k+1)  ) {
            amat[i][j][k] = 1;
        } else if (AMAT(i,j,k-1) > 1				   ||
                   V_CELL(0,i,j,k) != V_CELL(0,i,j,k-1)  ) {
            amat[i][j][k] = 1;
        }
    }
	free(v);
	for (i=0; i<8; i++) {
		v_cell[i]=free_3Dmatrix(v_cell[i]);
	}
	free(v_cell);
}
