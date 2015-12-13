#ifndef GEOMETRY_H
#define GEOMETRY_H

bool in_domain(float, float, float);
float *CrossProduct(float *, float *);
bool SameSide(float *,float *,float *,float *);
bool PointInTriangle(float *,float *,float *,float *);
double porosity(int ***);

#endif
