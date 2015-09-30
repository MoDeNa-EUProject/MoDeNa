/**
@cond

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2015 MoDeNa Consortium, All rights reserved.

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
    
@endcond
@file

MoDeNa low-level interface library

@author Henrik Rusche
@copyright  2014-2015, MoDeNa Project. GNU Public License.
*/

#ifndef __MODEL_H__
#define __MODEL_H__

#include "Python.h"
#include <stdbool.h>
#include "function.h"
#include "inputsoutputs.h"

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

extern PyTypeObject modena_model_tType;

extern PyObject *modena_SurrogateModel;

/**
@addtogroup C_interface_library
@{
*/

/**
stores a model and mapping for substitution
*/
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

/**
stores a surrogate model
*/
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

    struct modena_function_t *mf;

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

modena_model_t *modena_model_new
(
    const char *modelId
);

size_t modena_model_inputs_argPos
(
    const modena_model_t *self,
    const char *name
);

void modena_model_argPos_check(const modena_model_t *self);

size_t modena_model_inherited_inputs_argPos
(
    const modena_model_t *self,
    const char *name
);

size_t modena_model_outputs_argPos
(
    const modena_model_t *self,
    const char *name
);

size_t modena_model_inputs_size(const modena_model_t *self);

size_t modena_model_inherited_inputs_size(const modena_model_t *self);

size_t modena_model_outputs_size(const modena_model_t *self);

void modena_model_inputs_siunits
(
    const modena_model_t *self,
    const size_t i,
    modena_siunits_t *units
);

void modena_model_inherited_inputs_siunits
(
    const modena_model_t *self,
    const size_t i,
    modena_siunits_t *units
);

void modena_model_outputs_siunits
(
    const modena_model_t *self,
    const size_t i,
    modena_siunits_t *units
);

/**
modena_model_call returns:

- 201: requesting exit for new DOE without Restart
- 200: requesting exit for new DOE with Restart
- 100: updated model parameters, requesting to continue this run
- 1: failure
- 0: okay

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

/** @} */ // end of C_interface_library

__END_DECLS

#endif /* __MODEL_H__ */

