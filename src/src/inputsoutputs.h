/*
@cond

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2016 MoDeNa Consortium, All rights reserved.

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

Low-level interface library

@author Henrik Rusche
@copyright  2014-2016, MoDeNa Project. GNU Public License.
*/

#include "inline.h"
#include <stddef.h>

#ifndef __INPUTSOUTPUTS_H__
#define __INPUTSOUTPUTS_H__

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

/**
@addtogroup C_interface_library
@{
*/

// Forward declaration
struct modena_model_t;

/**
stores the exponents of the 7 SI units
*/
typedef struct
{
    int exponents[6];

} modena_siunits_t;

/**
stores the input data - both direct and inherited
*/
typedef struct
{
    double *inputs;

    double *inherited_inputs;

} modena_inputs_t;

/**
stores the outputs of a surrogate function
*/
typedef struct
{
    double *outputs;

} modena_outputs_t;

/**
creates from nothing
*/
modena_siunits_t *modena_siunits_new();

int modena_siunits_get(const modena_siunits_t *self, const size_t i);

void modena_siunits_destroy(modena_siunits_t *self);

modena_inputs_t *modena_inputs_new(const struct modena_model_t *self);

modena_outputs_t *modena_outputs_new(const struct modena_model_t *self);

void modena_inputs_destroy(modena_inputs_t *inputs);

void modena_outputs_destroy(modena_outputs_t *outputs);

INLINE_DECL void modena_inputs_set(modena_inputs_t *self, const size_t i, double x);

INLINE_DECL void modena_inherited_inputs_set
(
    modena_inputs_t *self,
    const size_t i,
    double x
);

INLINE_DECL double modena_inputs_get(const modena_inputs_t *self, const size_t i);

INLINE_DECL double modena_inherited_inputs_get
(
    const modena_inputs_t *self,
    const size_t i
);

INLINE_DECL double modena_outputs_get(const modena_outputs_t *self, const size_t i);

#ifdef HAVE_INLINE

INLINE_FUN void modena_inputs_set(modena_inputs_t *self, const size_t i, double x)
{
    self->inputs[i] = x;
}

INLINE_FUN void modena_inherited_inputs_set
(
    modena_inputs_t *self,
    const size_t i,
    double x
)
{
    self->inherited_inputs[i] = x;
}

INLINE_FUN double modena_inputs_get(const modena_inputs_t *self, const size_t i)
{
    return self->inputs[i];
}

INLINE_FUN double modena_inherited_inputs_get
(
    const modena_inputs_t *self,
    const size_t i
)
{
    return self->inherited_inputs[i];
}

INLINE_FUN double modena_outputs_get(const modena_outputs_t *self, const size_t i)
{
    return self->outputs[i];
}

#endif /* HAVE_INLINE */

/** @} */ // end of C_interface_library

__END_DECLS

#endif /* __INPUTSOUTPUTS_H__ */


