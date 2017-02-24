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

#ifndef __INDEXSET_H__
#define __INDEXSET_H__

#include "Python.h"

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

extern PyTypeObject modena_index_set_tType;

extern PyObject *modena_IndexSet;

/**
@addtogroup C_interface_library
@{
*/

/**
stores an index set
*/
typedef struct modena_index_set_t
{
    PyObject_HEAD

    PyObject *pIndexSet;

} modena_index_set_t;

modena_index_set_t *modena_index_set_new
(
    const char *indexSetId
);

size_t modena_index_set_get_index
(
    const modena_index_set_t *self,
    const char* name
);

const char* modena_index_set_get_name
(
    const modena_index_set_t *self,
    const size_t index
);

size_t modena_index_set_iterator_start
(
    const modena_index_set_t *self
);

size_t modena_index_set_iterator_end
(
    const modena_index_set_t *self
);

void modena_index_set_destroy(modena_index_set_t *indexSet);

/** @} */ // end of C_interface_library

__END_DECLS

#endif /* __INDEXSET_H__ */

