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

#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include "Python.h"
#include <ltdl.h>
#include "indexset.h"

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

// Forward declaration
struct modena_model_t;

extern PyTypeObject modena_function_tType;

extern PyObject *modena_SurrogateFunction;

/**
@addtogroup C_interface_library
@{
*/

/**
stores a surrogate function
*/
typedef struct modena_function_t
{
    PyObject_HEAD

    PyObject *pFunction;

    size_t parameters_size;

    lt_dlhandle handle;

    void (*function)
    (
        const struct modena_model_t* model,
        const double* i,
        double *o
    );

} modena_function_t;

modena_function_t *modena_function_new
(
    const char *functionId
);

modena_function_t *modena_function_new_from_model
(
    const struct modena_model_t *self
);

modena_index_set_t *modena_function_get_index_set
(
    const modena_function_t* self,
    const char* name
);

void modena_function_destroy(modena_function_t *model);

/** @} */ // end of C_interface_library

__END_DECLS

#endif /* __FUNCTION_H__ */

