/**
@file bubbleRadius.h
@brief functions related to the calculations of bubble radius.
@fn bool isNaN(double var)
@brief chech if the passed argument is not a number
@param double var - input variable
@return true, if the argument is NaN
@fn double bubbleRadius (const double m0, const double m1)
@brief radius of bubbles based on the moments 
@param const double m0 - moment of order zero
@param const double m1 - moment of order one
@return bubble radius
@fn double nodeRadius(const double &v)
@brief radius of bubbles at each node
@param const double v - volume of bubble
@return radius of bubble at the node
*/
bool isNaN(double var);
double bubbleRadius (const double m0, const double m1);
double nodeRadius(const double &v);

bool isNaN(double var)
{
    volatile double d = var;
    return d != d;
}
double bubbleRadius (const double m0, const double m1)
{
    double R;
    R   = pow((3.0*m1/(4.0*M_PI*m0)), 1.0/3.0);
    if (isNaN(R)) {
        R=init_size;
    }
    return R;
}

double nodeRadius(const double &v)
{
    return (pow((6.0*v/M_PI), 1/3.0));
}