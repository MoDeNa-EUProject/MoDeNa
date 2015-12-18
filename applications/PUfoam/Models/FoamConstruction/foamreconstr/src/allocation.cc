/*! \file
	\brief Shorthand functions for allocation of matrices
	\author Pavel Ferkl
	\author Juraj Kosek
	\ingroup foam_constr
*/
#include <stdlib.h>
//! allocate 3D double matrix
double ***alloc_3Ddmatrix (\
	int nx /**< [in] number of elements in 1st direction */, \
	int ny /**< [in] number of elements in 2nd direction */, \
	int nz /**< [in] number of elements in 3rd direction */)
{
	int i, j;
	double ***amat = NULL;

	amat = (double ***)calloc((size_t)nx, sizeof(double **));
	if (amat == NULL)
		return NULL;

	*amat = (double **)calloc((size_t)(nx*ny), sizeof(double *));
	if (*amat == NULL)
		return NULL;
	for (i = 1; i < nx; i++)
		amat[i] = &amat[0][i*ny];

	**amat = (double *)calloc((size_t)(nx*ny*nz), sizeof(double));
	if (**amat == NULL)
		return NULL;
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			amat[i][j] = &amat[0][0][i*ny*nz+j*nz];

	return amat;
}
//! free 3D double matrix
double ***free_3Ddmatrix (double ***amat /**< [in,out] matrix */)
{
	free ((void *)**amat);
	free ((void *) *amat);
	free ((void *)  amat);
	amat = NULL;
	return amat;
}
//! allocate 3D float matrix
float ***alloc_3Dfmatrix (\
	int nx /**< [in] number of elements in 1st direction */, \
	int ny /**< [in] number of elements in 2nd direction */, \
	int nz /**< [in] number of elements in 3rd direction */)
{
	int i, j;
	float ***amat = NULL;

	amat = (float ***)calloc((size_t)nx, sizeof(float **));
	if (amat == NULL)
		return NULL;

	*amat = (float **)calloc((size_t)(nx*ny), sizeof(float *));
	if (*amat == NULL)
		return NULL;
	for (i = 1; i < nx; i++)
		amat[i] = &amat[0][i*ny];

	**amat = (float *)calloc((size_t)(nx*ny*nz), sizeof(float));
	if (**amat == NULL)
		return NULL;
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			amat[i][j] = &amat[0][0][i*ny*nz+j*nz];

	return amat;
}
//! free 3D float matrix
float ***free_3Dfmatrix (float ***amat /**< [in,out] matrix */)
{
	free ((void *)**amat);
	free ((void *) *amat);
	free ((void *)  amat);
	amat = NULL;
	return amat;
}
//! allocate 3D integer matrix
int ***alloc_3Dmatrix (\
	int nx /**< [in] number of elements in 1st direction */, \
	int ny /**< [in] number of elements in 2nd direction */, \
	int nz /**< [in] number of elements in 3rd direction */)
{
	int i, j;
	int ***amat = NULL;

	amat = (int ***)calloc((size_t)nx, sizeof(int **));
	if (amat == NULL)
		return NULL;

	*amat = (int **)calloc((size_t)(nx*ny), sizeof(int *));
	if (*amat == NULL)
		return NULL;
	for (i = 1; i < nx; i++)
		amat[i] = &amat[0][i*ny];

	**amat = (int *)calloc((size_t)(nx*ny*nz), sizeof(int));
	if (**amat == NULL)
		return NULL;
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			amat[i][j] = &amat[0][0][i*ny*nz+j*nz];

	return amat;
}
//! free 3D integer matrix
int ***free_3Dmatrix (int ***amat /**< [in,out] matrix */)
{
	free ((void *)**amat);
	free ((void *) *amat);
	free ((void *)  amat);
	amat = NULL;
	return amat;
}
//! allocate 2D double matrix
double **alloc_dmatrix (\
	int nx /**< [in] number of elements in 1st direction */, \
	int ny /**< [in] number of elements in 2nd direction */)
{
	int i;
	double **amat = NULL;

	/*
	 * Allocate the matrix.
	 */
	amat = (double **)calloc((size_t)nx, sizeof(double *));
	if (amat == NULL)
		return NULL;

	*amat = (double *)calloc((size_t)(nx*ny), sizeof(double));
	if (*amat == NULL)
		return NULL;

	for (i = 1; i < nx; i++)
		amat[i] = &amat[0][i*ny];

	return amat;
}
//! free 2D double matrix
double **free_dmatrix (double **amat /**< [in,out] matrix */)
{
	free ((void *)*amat);
	free ((void *) amat);
	amat = NULL;
	return amat;
}
//! allocate 2D integer matrix
float **alloc_fmatrix (\
	int nx /**< [in] number of elements in 1st direction */, \
	int ny /**< [in] number of elements in 2nd direction */)
{
	int i;
	float **amat = NULL;

	/*
	 * Allocate the matrix.
	 */
	amat = (float **)calloc((size_t)nx, sizeof(float *));
	if (amat == NULL)
		return NULL;

	*amat = (float *)calloc((size_t)(nx*ny), sizeof(float));
	if (*amat == NULL)
		return NULL;

	for (i = 1; i < nx; i++)
		amat[i] = &amat[0][i*ny];

	return amat;
}
//! free 2D integer matrix
float **free_fmatrix (float **amat /**< [in,out] matrix */)
{
	free ((void *)*amat);
	free ((void *) amat);
	amat = NULL;
	return amat;
}

int **alloc_matrix (\
	int nx /**< [in] number of elements in 1st direction */, \
	int ny /**< [in] number of elements in 2nd direction */)
{
	int i;
	int **amat = NULL;

	/*
	 * Allocate the matrix.
	 */
	amat = (int **)calloc((size_t)nx, sizeof(int *));
	if (amat == NULL)
		return NULL;

	*amat = (int *)calloc((size_t)(nx*ny), sizeof(int));
	if (*amat == NULL)
		return NULL;

	for (i = 1; i < nx; i++)
		amat[i] = &amat[0][i*ny];

	return amat;
}

int **free_matrix (int **amat /**< [in,out] matrix */)
{
	free ((void *)*amat);
	free ((void *) amat);
	amat = NULL;
	return amat;
}
