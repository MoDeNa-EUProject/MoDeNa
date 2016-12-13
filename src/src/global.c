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

#undef __INLINE_H__
#define __INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif

#include "global.h"

#ifdef HAVE_INLINE       /* disable compilation of definitions in .h files */
#undef HAVE_INLINE
#endif

#include "indexset.h"
#include "function.h"
#include "model.h"
#include <execinfo.h>

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


// Initialise global variable
thread_local int modena_error_code = 0;

PyObject *modena_DoesNotExist = NULL;
PyObject *modena_OutOfBounds = NULL;
PyObject *modena_ParametersNotValid = NULL;

struct modena_errordesc
{
    int  code;
    const char *message;
} modena_errordesc[] =
{
    { MODENA_SUCCESS, "No error" },
    { MODENA_MODEL_NOT_FOUND, "Surrogate model not found in database" },
    { MODENA_FUNCTION_NOT_FOUND, "Surrogate function not found in database" },
    { MODENA_INDEX_SET_NOT_FOUND, "Index set not found in database" }
};

const char* modena_error_message(int error_code)
{
    return modena_errordesc[modena_error_code].message;
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC    /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

// TODO: Support non-Gcc compilers here
PyMODINIT_FUNC initlibmodena(void) __attribute__((constructor));

PyMODINIT_FUNC initlibmodena(void)
{
    // Initialize the Python Interpreter
    if(!Py_IsInitialized())
    {
        Py_Initialize();
    }

    if(PyType_Ready(&modena_index_set_tType) < 0)
    {
        return;
    }

    if(PyType_Ready(&modena_function_tType) < 0)
    {
        return;
    }

    if(PyType_Ready(&modena_model_tType) < 0)
    {
        return;
    }

    PyObject* module = Py_InitModule3
    (
        "libmodena",
        module_methods,
        "Module that creates extension types for modena framework."
    );

    if(!module)
    {
        return;
    }

    if(!modena_DoesNotExist)
    {
        Py_INCREF(&modena_index_set_tType);
        PyModule_AddObject
        (
            module,
            "modena_index_set_t",
            (PyObject *) &modena_index_set_tType
        );

        Py_INCREF(&modena_function_tType);
        PyModule_AddObject
        (
            module,
            "modena_function_t",
            (PyObject *) &modena_function_tType
        );

        Py_INCREF(&modena_model_tType);
        PyModule_AddObject
        (
            module,
            "modena_model_t",
            (PyObject *) &modena_model_tType
        );

        PyObject *pName = PyString_FromString("modena.SurrogateModel");
        if(!pName){ Modena_PyErr_Print(); }

        PyObject *pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if(!pModule){ Modena_PyErr_Print(); }

        PyObject *pDict = PyModule_GetDict(pModule); // Borrowed ref
        if(!pDict){ Modena_PyErr_Print(); }


        pName = PyString_FromString("IndexSet");
        if(!pName){ Modena_PyErr_Print(); }

        modena_IndexSet = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!modena_IndexSet){ Modena_PyErr_Print(); }

        pName = PyString_FromString("SurrogateFunction");
        if(!pName){ Modena_PyErr_Print(); }

        modena_SurrogateFunction = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!modena_SurrogateFunction){ Modena_PyErr_Print(); }

        pName = PyString_FromString("SurrogateModel");
        if(!pName){ Modena_PyErr_Print(); }

        modena_SurrogateModel = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!modena_SurrogateModel){ Modena_PyErr_Print(); }


        pName = PyString_FromString("DoesNotExist");
        if(!pName){ Modena_PyErr_Print(); }

        modena_DoesNotExist = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!modena_DoesNotExist){ Modena_PyErr_Print(); }

        pName = PyString_FromString("ParametersNotValid");
        if(!pName){ Modena_PyErr_Print(); }

        modena_ParametersNotValid = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!modena_ParametersNotValid){ Modena_PyErr_Print(); }

        pName = PyString_FromString("OutOfBounds");
        if(!pName){ Modena_PyErr_Print(); }

        modena_OutOfBounds = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!modena_OutOfBounds){ Modena_PyErr_Print(); }


        Py_DECREF(pModule);
    }
}


void modena_print_backtrace()
{
    void* tracePtrs[100];
    int count = backtrace( tracePtrs, 100 );

    char** funcNames = backtrace_symbols( tracePtrs, count );
    // Print the stack trace
    int ii;
    for( ii = 0; ii < count; ii++ )
        printf( "%s\n", funcNames[ii] );

    // Free the string pointers
    free( funcNames );

    exit(1);
}

