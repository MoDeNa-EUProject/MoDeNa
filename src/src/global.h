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

Description
    Interface Library

Authors
    Henrik Rusche

Contributors
*/

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include "Python.h"
#include "inline.h"
#include <stdbool.h>

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

#ifndef thread_local
# if __STDC_VERSION__ >= 201112 && !defined __STDC_NO_THREADS__
#  define thread_local _Thread_local
# elif defined _WIN32 && ( \
       defined _MSC_VER || \
       defined __ICL || \
       defined __DMC__ || \
       defined __BORLANDC__ )
#  define thread_local __declspec(thread)
/* note that ICC (linux) and Clang are covered by __GNUC__ */
# elif defined __GNUC__ || \
       defined __SUNPRO_C || \
       defined __xlC__
#  define thread_local __thread
# else
#  error "Cannot define thread_local"
# endif
#endif

// Declare variable for error handling
extern thread_local int modena_error_code;

// Declare python exceptions
extern PyObject *modena_DoesNotExist;
extern PyObject *modena_OutOfBounds;
extern PyObject *modena_ParametersNotValid;

enum modena_error_t
{
    MODENA_SUCCESS,
    MODENA_MODEL_NOT_FOUND,
    MODENA_FUNCTION_NOT_FOUND,
    MODENA_INDEX_SET_NOT_FOUND,
    MODENA_MODEL_LAST
};

// Returns true when an error has been raised
INLINE_DECL bool modena_error_occurred();

// Returns error code and resets it
INLINE_DECL int modena_error();

#ifdef HAVE_INLINE

INLINE_FUN bool modena_error_occurred()
{
    return modena_error_code != MODENA_SUCCESS;
}

INLINE_FUN int modena_error()
{
    int ret = modena_error_code;
    modena_error_code = 0;
    return ret;
}

#endif /* HAVE_INLINE */

// Returns error message for error code
const char* modena_error_message(int error_code);

void modena_print_backtrace();

#define Modena_Info_Print(...)                                                \
    char Modena_message[256];                                                 \
    sprintf(Modena_message, __VA_ARGS__);                                     \
    fprintf(stdout, "%s in line %i of %s\n", Modena_message,  __LINE__, __FILE__);

#define Modena_Error_Print(...)                                               \
    char Modena_message[256];                                                 \
    sprintf(Modena_message, __VA_ARGS__);                                     \
    fprintf(stderr, "%s in line %i of %s\n", Modena_message, __LINE__, __FILE__);

#define Modena_PyErr_Print()                                                  \
    PyErr_Print();                                                            \
    Modena_Error_Print("Error in python catched");                            \
    modena_print_backtrace();

__END_DECLS

#endif /* __GLOBAL_H__ */

