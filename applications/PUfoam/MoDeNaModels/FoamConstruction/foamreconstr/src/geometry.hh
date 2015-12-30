#ifndef GEOMETRY_H
#define GEOMETRY_H

bool in_domain(double, double, double);
double *CrossProduct(double *, double *);
bool SameSide(double *,double *,double *,double *);
bool PointInTriangle(double *,double *,double *,double *);
double porosity(int ***);

#endif
