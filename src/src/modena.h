/*

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014 MoDeNa Consortium, All rights reserved.

License
    This file is part of Modena.

    The Modena interface library is free software; you can redistribute it
    and/or modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

Description
    Interface Library

Authors
    Henrik Rusche

Contributors
*/

#ifndef __MODENA_H__
#define __MODENA_H__

#include "Python.h"
#include <stdbool.h>
#include <ltdl.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum modena_error_t
{
    MODENA_SUCCESS,
    MODENA_MODEL_NOT_FOUND,
    MODENA_MODEL_LAST
};

struct modena_errordesc
{
    int  code;
    const char *message;
} modena_errordesc[] =
{
    { MODENA_SUCCESS, "No error" },
    { MODENA_MODEL_NOT_FOUND, "Surrogate model not found in database" }
};

// Returns error and resets it
int modena_error();

// Returns true when an error has been raised
bool modena_error_occurred();

// Returns error message for error code
const char* modena_error_message(int error_code);

#define Modena_PyErr_Print() \
    fprintf(stderr, "Error in python catched in %i of %s\n", __LINE__, __FILE__); \
    PyErr_Print(); \
    exit(1);

// modena_units_t stores the exponents of the 7 SI units
typedef struct
{
    int exponents[6];

} modena_siunits_t;

// modena_inputs_t stores the input data - both direct and inherited
typedef struct
{
	double *inputs;

	double *inherited_inputs;

} modena_inputs_t;

// modena_units_t stores the outputs of a surrogate function
typedef struct
{
	double *outputs;

} modena_outputs_t;

// modena_function_t stores a surrogate function
typedef struct
{
    lt_dlhandle handle;

    void (*function)
    (
        const double* p,
        const double* in_i,
        const double* i,
        double *o
    );

} modena_function_t;

// modena_substitute_model_t stores a model and mapping for substitution
typedef struct modena_substitute_model_t
{
    struct modena_model_t *model;

    modena_inputs_t *inputs;

    modena_outputs_t *outputs;

    size_t map_inputs_size;

    size_t *map_inputs;

    size_t map_outputs_size;

    size_t *map_outputs;

} modena_substitute_model_t;

// modena_function_t stores a surrogate model
typedef struct modena_model_t
{
    PyObject_HEAD;

    PyObject *pModel;

    size_t outputs_size;

    size_t inputs_size;

    size_t inputs_minMax_size;

    double *inputs_min;

    double *inputs_max;

    bool *argPos_used;

    size_t inherited_inputs_size;

    size_t parameters_size;

    double *parameters;

    modena_function_t *mf;

    void (*function)
    (
        const double* p,
        const double* in_i,
        const double* i,
        double *o
    );

    size_t substituteModels_size;

    modena_substitute_model_t *substituteModels;

} modena_model_t;

modena_siunits_t *modena_siunits_new();

int modena_siunits_get(const modena_siunits_t *units, const size_t i)
{
    return units->exponents[i];
}

void modena_siunits_destroy(modena_siunits_t *units);

modena_inputs_t *modena_inputs_new(const modena_model_t *model);

modena_outputs_t *modena_outputs_new(const modena_model_t *model);

void modena_inputs_set(modena_inputs_t *inputs, const size_t i, double x)
{
    inputs->inputs[i] = x;
}

void modena_inherited_inputs_set
(
    modena_inputs_t *inputs,
    const size_t i,
    double x
)
{
    inputs->inherited_inputs[i] = x;
}

double modena_inputs_get(const modena_inputs_t *inputs, const size_t i)
{
    return inputs->inputs[i];
}

double modena_inherited_inputs_get(const modena_inputs_t *inputs, const size_t i)
{
    return inputs->inherited_inputs[i];
}

double modena_outputs_get(const modena_outputs_t *outputs, const size_t i)
{
    return outputs->outputs[i];
}

void modena_inputs_destroy(modena_inputs_t *inputs);

void modena_outputs_destroy(modena_outputs_t *outputs);

modena_model_t *modena_model_new
(
    const char *modelId
);

size_t modena_model_inputs_argPos(const modena_model_t *model, const char *name);

void modena_model_argPos_check(const modena_model_t *model);

size_t modena_model_inherited_inputs_argPos
(
    const modena_model_t *model,
    const char *name
);

size_t modena_model_outputs_argPos(const modena_model_t *model, const char *name);

size_t modena_model_inputs_size(const modena_model_t *model);

size_t modena_model_inherited_inputs_size(const modena_model_t *model);

size_t modena_model_outputs_size(const modena_model_t *model);

void modena_model_inputs_siunits
(
    const modena_model_t *model,
    const size_t i,
    modena_siunits_t *units
);

void modena_model_inherited_inputs_siunits
(
    const modena_model_t *model,
    const size_t i,
    modena_siunits_t *units
);

void modena_model_outputs_siunits
(
    const modena_model_t *model,
    const size_t i,
    modena_siunits_t *units
);

/*
modena_model_call returns:

201: requesting exit for new DOE without Restart
200: requesting exit for new DOE with Restart
100: updated model parameters, requesting to continue this run
1: failure
0: okay

If exit is requested, do what's necessary and exit with the same error code!
*/
int modena_model_call
(
    modena_model_t *model,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
);

void modena_model_call_no_check
(
    modena_model_t *model,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
);

void modena_model_destroy(modena_model_t *model);

modena_function_t *modena_function_new
(
    const modena_model_t *m
);

void modena_function_destroy(modena_function_t *model);


__END_DECLS

#endif /* __MODENA_H__ */

