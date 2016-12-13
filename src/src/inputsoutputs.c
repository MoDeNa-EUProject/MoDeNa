/*

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
*/

#include "global.h"

#undef __INLINE_H__
#define __INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif

#include "inputsoutputs.h"

#ifdef HAVE_INLINE       /* disable compilation of definitions in .h files */
#undef HAVE_INLINE
#endif

#include "model.h"

modena_siunits_t *modena_siunits_new()
{
    return malloc(sizeof(modena_siunits_t));
}

void modena_siunits_destroy(modena_siunits_t *self)
{
    free(self);
}

modena_inputs_t *modena_inputs_new(const modena_model_t *self)
{
    modena_inputs_t *i = malloc(sizeof(modena_inputs_t));
    i->inputs = malloc(self->inputs_size*sizeof(double));
    return i;
}

modena_outputs_t *modena_outputs_new(const modena_model_t *self)
{
    modena_outputs_t *o = malloc(sizeof(modena_outputs_t));
    o->outputs = malloc(self->outputs_size*sizeof(double));
    return o;
}

void modena_inputs_destroy(modena_inputs_t *self)
{
    free(self->inputs);
    free(self);
}

void modena_outputs_destroy(modena_outputs_t *self)
{
    free(self->outputs);
    free(self);
}

