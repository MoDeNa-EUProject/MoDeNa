#ifndef ALLOCATION_H
#define ALLOCATION_H

float ***alloc_3Dfmatrix (int, int, int);
float ***free_3Dfmatrix  (float ***);
int ***alloc_3Dmatrix (int, int, int);
int ***free_3Dmatrix  (int ***);
double **alloc_dmatrix (int, int);
double **free_dmatrix  (double **);
float **alloc_fmatrix (int, int);
float **free_fmatrix  (float **);
int    **alloc_matrix  (int, int);
int    **free_matrix   (int **);

#endif
