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

#include "function.h"
#include "structmember.h"
#include "global.h"
#include "model.h"

PyObject *modena_SurrogateFunction = NULL;

void modena_function_load_library(modena_function_t* self)
{
    PyObject *pFunctionName =
        PyObject_GetAttrString(self->pFunction, "functionName");
    if(!pFunctionName){ Modena_PyErr_Print(); }

    PyObject *pLibraryName =
        PyObject_GetAttrString(self->pFunction, "libraryName");
    if(!pLibraryName){ Modena_PyErr_Print(); }

    self->handle = lt_dlopen(PyString_AsString(pLibraryName));

    if(!self->handle)
    {
        Modena_Error_Print
        (
           "lt_dlopen: Could not open library %s\nlt_dlopen: %s",
           PyString_AsString(pLibraryName),
           lt_dlerror()
        );
        exit(1);
    }

    self->function = lt_dlsym(self->handle, PyString_AsString(pFunctionName));
    if(!self->function)
    {
        Modena_Error_Print
        (
            "lt_dlsym: Could not find function %s in library %s"
            "lt_dlsym: %s",
            PyString_AsString(pFunctionName),
            PyString_AsString(pLibraryName),
            lt_dlerror()
        );
        lt_dlclose(self->handle);
        exit(1);
    }

    Py_DECREF(pFunctionName);
    Py_DECREF(pLibraryName);

    PyObject *pParameters =
        PyObject_GetAttrString(self->pFunction, "parameters");
    if(!pParameters){ Modena_PyErr_Print(); }
    self->parameters_size = PyObject_Size(pParameters);
    Py_DECREF(pParameters);
}


modena_function_t *modena_function_new
(
    const char *functionId
)
{
    // Initialize the Python Interpreter
    if(!Py_IsInitialized())
    {
        Py_Initialize();
    }

    // Initialize this module
    initlibmodena();

    PyObject *args = PyTuple_New(0);
    PyObject *kw = Py_BuildValue("{s:s}", "functionId", functionId);

    PyObject *pNewObj = PyObject_Call
    (
        (PyObject *) &modena_function_tType,
        args,
        kw
    );

    Py_DECREF(args);
    Py_DECREF(kw);
    if(!pNewObj)
    {
        if(PyErr_ExceptionMatches(modena_DoesNotExist))
        {
            PyErr_Clear();

            PyObject *pRet = PyObject_CallMethod
            (
                modena_SurrogateFunction,
                "exceptionLoad",
                "(z)",
                functionId
            );
            if(!pRet){ Modena_PyErr_Print(); }
            int ret = PyInt_AsLong(pRet);
            Py_DECREF(pRet);

            modena_error_code = ret;
            return NULL;
        }
        else
        {
            Modena_PyErr_Print();
        }
    }


    modena_function_t *self = (modena_function_t *) pNewObj;

    if(lt_dlinit())
    {
        Modena_Error_Print("lt_dlinit: %s", lt_dlerror());
        exit(1);
    }

    modena_function_load_library(self);

    return self;
}

modena_function_t *modena_function_new_from_model
(
    const modena_model_t *m
)
{
    modena_function_t *self = malloc(sizeof(modena_function_t));

    if(lt_dlinit())
    {
        Modena_Error_Print("lt_dlinit: %s", lt_dlerror());
        exit(1);
    }

    self->pFunction =
        PyObject_GetAttrString(m->pModel, "surrogateFunction");
    if(!self->pFunction){ Modena_PyErr_Print(); }

    modena_function_load_library(self);

    return self;
}

modena_index_set_t *modena_function_get_index_set
(
    const modena_function_t* self,
    const char* name
)
{
    PyObject *pRet = PyObject_CallMethod
    (
        self->pFunction,
        "indexSet",
        "(z)",
        name
    );
    if(!pRet){ Modena_PyErr_Print(); }

    PyObject *args = PyTuple_New(0);
    PyObject *kw = Py_BuildValue("{s:O}", "indexSet", pRet);

    PyObject *pNewObj = PyObject_Call
    (
        (PyObject *) &modena_index_set_tType,
        args,
        kw
    );

    Py_DECREF(args);
    Py_DECREF(kw);
    if(!pNewObj)
    {
       Modena_PyErr_Print();
    }

    return (modena_index_set_t *) pNewObj;
}

/* Destructor, deallocates the memory block occupied by a surrogate function
 */
void modena_function_destroy(modena_function_t *self)
{
    lt_dlclose(self->handle);
    lt_dlexit();
    Py_XDECREF(self->pFunction);
    free(self);
}

/* C-Python: Destructor, exposed as the __del__ method in Python.
 */
static void modena_function_t_dealloc(modena_function_t* self)
{
    modena_function_destroy(self);
}

/* C-Python: Member-Table
 *
 * Structure which describes an attribute of a type which corresponds to a C
 * struct member. Its fields are:
 *
 * Field  C Type       Meaning
 * ------ ----------  --------------------------------------------------------
 * name   char *      name of the member
 * type   int         the type of the member in the C struct
 * offset Py_ssize_t  the offset in bytes that the member is located on the
 *                    type's object struct
 * flags  int         flag bits indicating if the field should be read-only or 
 *                    writable
 * doc    char *      points to the contents of the docstring
 */
static PyMemberDef modena_function_t_members[] = {
    {NULL}  /* Sentinel */
};

/* C-Python: Method-Table
 *
 * Structure used to describe a method of an extension type. This structure has
 * four fields:
 *
 * Field     C Type       Meaning
 * -------   -----------  ----------------------------------------------------
 * ml_name   char *       name of the method
 * ml_meth   PyCFunction  pointer to the C implementation
 * ml_flags  int          flag bits indicating how the call should be
 *                        constructed
 * ml_doc    char *       points to the contents of the docstring
 */
static PyMethodDef modena_function_t_methods[] = {
    {NULL}  /* Sentinel */
};

/* C-Python: Initialiser, exposed in Python as the method: __init__
 */
static int modena_function_t_init
(
   modena_function_t *self,
   PyObject *args,
   PyObject *kwds
)
{
    PyObject *pFunction=NULL;
    char *functionId=NULL;

    static char *kwlist[] = {"function", "functionId", NULL};

    if
    (
        !PyArg_ParseTupleAndKeywords
        (
            args,
            kwds,
            "|Os",
            kwlist,
            &pFunction,
            &functionId
        )
    )
    {
        Modena_PyErr_Print();
    }

    if(!pFunction)
    {
        self->pFunction = PyObject_CallMethod
        (
            modena_SurrogateFunction,
            "load",
            "(z)",
            functionId
        );

        if(!self->pFunction)
        {
            PyErr_SetString(modena_DoesNotExist, "Function does not exist");

            Modena_PyErr_Print();
        }
    }
    else
    {
        Py_INCREF(pFunction);
        self->pFunction = pFunction;
    }

    // Shouldn't I load the function from the library here?

    PyObject *pParameters =
        PyObject_GetAttrString(self->pFunction, "parameters");
    if(!pParameters){ Modena_PyErr_Print(); }
    self->parameters_size = PyObject_Size(pParameters);
    Py_DECREF(pParameters);

    return 0;
}

/* C-Python: Constructor, exposed in Python as the method: __new__
 */
static PyObject *modena_function_t_new
(
    PyTypeObject *type,
    PyObject *args,
    PyObject *kwds
)
{
    modena_function_t *self;

    self = (modena_function_t *)type->tp_alloc(type, 0);
    if(self)
    {
        self->pFunction = NULL;
        self->parameters_size = 0;
    }

    return (PyObject *)self;
}

/* C-Python: The C structure used to describe the modena_model type.
 */
PyTypeObject modena_function_tType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "modena.modena_function_t", /*tp_name*/
    sizeof(modena_function_t), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)modena_function_t_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "modena_function_t objects", /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    modena_function_t_methods, /* tp_methods */
    modena_function_t_members, /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)modena_function_t_init, /* tp_init */
    0,                         /* tp_alloc */
    modena_function_t_new,  /* tp_new */
};
