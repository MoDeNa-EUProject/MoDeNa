#include <stdlib.h>
float ***alloc_3Dfmatrix (int nx, int ny, int nz)
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

float ***free_3Dfmatrix (float ***amat)
{
	free ((void *)**amat);
	free ((void *) *amat);
	free ((void *)  amat);
	amat = NULL;
	return amat;
}

int ***alloc_3Dmatrix (int nx, int ny, int nz)
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

int ***free_3Dmatrix (int ***amat)
{
	free ((void *)**amat);
	free ((void *) *amat);
	free ((void *)  amat);
	amat = NULL;
	return amat;
}

double **alloc_dmatrix (int nx, int ny)
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

double **free_dmatrix (double **amat)
{
	free ((void *)*amat);
	free ((void *) amat);
	amat = NULL;
	return amat;
}

float **alloc_fmatrix (int nx, int ny)
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

float **free_fmatrix (float **amat)
{
	free ((void *)*amat);
	free ((void *) amat);
	amat = NULL;
	return amat;
}

int **alloc_matrix (int nx, int ny)
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

int **free_matrix (int **amat)
{
	free ((void *)*amat);
	free ((void *) amat);
	amat = NULL;
	return amat;
}
