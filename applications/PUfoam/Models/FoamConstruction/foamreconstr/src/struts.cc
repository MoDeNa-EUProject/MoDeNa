/*! \file
	\brief Functions governing creation of struts
	\author Pavel Ferkl
	\ingroup foam_constr
*/
#include "globals.hh"
#include "geometry.hh"
#include "edges.hh"
#include "nodes.hh"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <iostream>
using namespace std;
using namespace globals;
//! Objective function for optimization. We want to have certain strut porosity.
double fn1 (\
	double x /**< [in] independent variable (`dstrut`) */,\
	void * p /**< [in] function parameters (`fn1_params`) */)
{
	// (void)(p); /* avoid unused parameter warning */
	struct fn1_params * params = (struct fn1_params *)p;
	int sv = (params->sv);
    int incmax = (params->incmax);
    int vmax = (params->vmax);
	double **vert = (params->vert);
	int **vinc = (params->vinc);
	int ***smat = (params->smat);
    bool report = (params->report);
	int i,j,k;
	dedge=x;
	for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	for (k = 0; k < nz; k++)
		smat[i][j][k] = 0;
	if (createEdges) {
		makeEdgeStruts(smat,sv,vert,vinc,report);
	}
	if (createNodes) {
		makeNodeStruts(smat,sv,incmax,vmax,vert,vinc,report);
	}
	double por = porosity(smat);
	return pow(por-strutPorosity,2);
}
//! GSL minimizer for univariate scalar function.
int optim (\
	void * params /**< [in] function parameters */,\
	double &m /**< [in] initial guess of minimum location */,\
	double &a /**< [in] lower bracket bound */,\
	double &b /**< [in] upper bracket bound */,\
	bool report /**< [in] show output */)
{
	int status;
	int iter = 0, max_iter = 100;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	gsl_function F;
    gsl_set_error_handler_off ();

	F.function = &fn1;
	F.params = params;

	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc (T);
	status = gsl_min_fminimizer_set (s, &F, m, a, b);
    if (status == GSL_EINVAL) {
        if (report) {
            cout << "endpoints do not enclose a minimum" << endl;
        }
        return(1);
    } else if (status != GSL_SUCCESS) {
        cout << "something is wrong in optim function" << endl;
    }

    if (report) {
    	printf ("using %s method\n", gsl_min_fminimizer_name (s));
    	printf ("%5s [%9s, %9s] %9s %9s\n",
            "iter", "lower", "upper", "min", "err(est)");
    	printf ("%5d [%.7f, %.7f] %.7f %.7f\n",iter, a, b, m, b - a);
    }

	do {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status
        = gsl_min_test_interval (a, b, 0.01, 0.0);

        if (status == GSL_SUCCESS)
        if (report) {
            printf ("Converged:\n");
            printf ("%5d [%.7f, %.7f] %.7f %.7f\n", iter, a, b, m, b - a);
        }
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_min_fminimizer_free (s);

	return status;
}
