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
*/

#include "modena.h"
#include "structmember.h"

// Forward declarations
static PyTypeObject modena_model_tType;

modena_siunits_t *modena_siunits_new()
{
    return malloc(sizeof(modena_siunits_t));
}

void modena_siunits_destroy(modena_siunits_t *units)
{
    free(units);
}

modena_inputs_t *modena_inputs_new(const modena_model_t *m)
{
    modena_inputs_t *i = malloc(sizeof(modena_inputs_t));
    i->inputs = malloc(m->inputs_size*sizeof(double));
    i->inherited_inputs = malloc
    (
        m->inherited_inputs_size*sizeof(double)
    );
    return i;
}

modena_outputs_t *modena_outputs_new(const modena_model_t *m)
{
    modena_outputs_t *o = malloc(sizeof(modena_outputs_t));
    o->outputs = malloc(m->outputs_size*sizeof(double));
    return o;
}

void modena_inputs_destroy(modena_inputs_t *i)
{
    free(i->inputs);
    free(i->inherited_inputs);
    free(i);
}

void modena_outputs_destroy(modena_outputs_t *o)
{
    free(o->outputs);
    free(o);
}

void modena_substitute_model_calculate_maps
(
    modena_substitute_model_t *sm,
    modena_model_t *parent
)
{
    PyObject *pMaps = PyObject_CallMethod
    (
        parent->pModel, "calculate_maps", "(O)", sm->model->pModel
    );
    if(!pMaps){ Modena_PyErr_Print(); }

    PyObject *pMapOutputs = PyTuple_GET_ITEM(pMaps, 0); // Borrowed ref
    if(!pMapOutputs){ Modena_PyErr_Print(); }
    PyObject *pSeq = PySequence_Fast(pMapOutputs, "expected a sequence");
    sm->map_outputs_size = PySequence_Size(pMapOutputs);
    sm->map_outputs = malloc(sm->map_outputs_size*sizeof(double));
    size_t i;
    for(i = 0; i < sm->map_outputs_size; i++)
    {
        sm->map_outputs[i] = PyInt_AsSsize_t(PyList_GET_ITEM(pSeq, i));
        parent->argPos_used[sm->map_outputs[i]] = true;
    }
    sm->map_outputs_size /= 2;
    Py_DECREF(pSeq);
    Py_DECREF(pMapOutputs);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    PyObject *pMapInputs = PyTuple_GET_ITEM(pMaps, 1); // Borrowed ref
    if(!pMapInputs){ Modena_PyErr_Print(); }
    pSeq = PySequence_Fast(pMapInputs, "expected a sequence");
    sm->map_inputs_size = PySequence_Size(pMapInputs);
    sm->map_inputs = malloc(sm->map_inputs_size*sizeof(double));
    for(i = 0; i < sm->map_inputs_size; i++)
    {
        sm->map_inputs[i] = PyInt_AsSsize_t(PyList_GET_ITEM(pSeq, i));
    }
    sm->map_inputs_size /= 2;
    Py_DECREF(pSeq);
    Py_DECREF(pMapInputs);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    Py_DECREF(pMaps);
}

void modena_model_read_substituteModels(modena_model_t *m)
{
    PyObject *pSubstituteModels = PyObject_GetAttrString
    (
        m->pModel, "substituteModels"
    );
    if(!pSubstituteModels){ Modena_PyErr_Print(); }

    PyObject *pSeq = PySequence_Fast
    (
        pSubstituteModels, "expected a sequence"
    );
    m->substituteModels_size = PySequence_Size(pSubstituteModels);
    m->substituteModels =
        malloc(m->substituteModels_size*sizeof(modena_substitute_model_t));
    size_t i;
    for(i = 0; i < m->substituteModels_size; i++)
    {
        PyObject *args = PyTuple_New(0);
        PyObject *kw = Py_BuildValue
        (
            "{s:O}", "model", PyList_GET_ITEM(pSeq, i)
        );

        m->substituteModels[i].model = (modena_model_t *) PyObject_Call
        (
            (PyObject *) &modena_model_tType, args, kw
        );
        Py_DECREF(args);
        Py_DECREF(kw);
        if(!m->substituteModels[i].model){ Modena_PyErr_Print(); }

        m->substituteModels[i].inputs = modena_inputs_new
        (
            m->substituteModels[i].model
        );

        m->substituteModels[i].outputs = modena_outputs_new
        (
            m->substituteModels[i].model
        );

        modena_substitute_model_calculate_maps(&m->substituteModels[i], m);
    }

    Py_DECREF(pSeq);
    Py_DECREF(pSubstituteModels);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }
}

void modena_model_get_minMax
(
    modena_model_t *m
)
{
    PyObject *pObj = PyObject_CallMethod(m->pModel, "minMax", NULL);
    if(!pObj){ Modena_PyErr_Print(); }

    PyObject *pMin = PyTuple_GET_ITEM(pObj, 0); // Borrowed ref
    PyObject *pSeq = PySequence_Fast(pMin, "expected a sequence");
    m->inputs_minMax_size = PySequence_Size(pMin);
    m->inputs_min = malloc(m->inputs_minMax_size*sizeof(double));
    size_t i;
    for(i = 0; i < m->inputs_minMax_size; i++)
    {
        m->inputs_min[i] = PyFloat_AsDouble(PyList_GET_ITEM(pSeq, i));
    }
    Py_DECREF(pSeq);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    PyObject *pMax = PyTuple_GET_ITEM(pObj, 1); // Borrowed ref
    pSeq = PySequence_Fast(pMax, "expected a sequence");
    m->inputs_max = malloc(m->inputs_minMax_size*sizeof(double));
    for(i = 0; i < m->inputs_minMax_size; i++)
    {
        m->inputs_max[i] = PyFloat_AsDouble(PyList_GET_ITEM(pSeq, i));
    }
    Py_DECREF(pSeq);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    Py_DECREF(pObj);
}

static void
modena_model_t_dealloc(modena_model_t* self)
{
    modena_model_destroy(self);
}

static PyObject *
modena_model_t_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    modena_model_t *self;

    self = (modena_model_t *)type->tp_alloc(type, 0);
    if (self != NULL)
    {
        // Init to null required?
    }

    return (PyObject *)self;
}

static int
modena_model_t_init(modena_model_t *self, PyObject *args, PyObject *kwds)
{
    PyObject *pParameters=NULL, *pModel=NULL;
    char *modelId=NULL;

    static char *kwlist[] = {"model", "modelId", "parameters"};

    if
    (
        !PyArg_ParseTupleAndKeywords
        (
            args, kwds, "|OsO", kwlist, &pModel, &modelId, &pParameters
        )
    )
    {
        return -1;
    }

    if(!pModel)
    {
        PyObject *pName = PyString_FromString("modena.SurrogateModel");
        if(!pName){ Modena_PyErr_Print(); }

        PyObject *pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if(!pModule){ Modena_PyErr_Print(); }

        PyObject *pDict = PyModule_GetDict(pModule); // Borrowed ref
        if(!pDict){ Modena_PyErr_Print(); }

        pName = PyString_FromString("SurrogateModel");
        if(!pName){ Modena_PyErr_Print(); }

        PyObject *sModel = PyObject_GetItem(pDict, pName);
        Py_DECREF(pName);
        if(!sModel){ Modena_PyErr_Print(); }

        self->pModel = PyObject_CallMethod(sModel, "load", "(z)", modelId);
        Py_DECREF(sModel);
        if(!self->pModel){ Modena_PyErr_Print(); }

        modena_model_get_minMax(self);

        Py_DECREF(pModule);
    }
    else
    {
        Py_INCREF(pModel);
        self->pModel = pModel;

        // Set everything to zero
        self->inputs_minMax_size = 0;
        self->inputs_min = NULL;
        self->inputs_max = NULL;
    }

    //PyObject_Print(self->pModel, stdout, 0);
    //printf("\n");

    PyObject *pOutputs = PyObject_GetAttrString(self->pModel, "outputs");
    if(!pOutputs){ Modena_PyErr_Print(); }
    self->outputs_size = PyDict_Size(pOutputs);
    Py_DECREF(pOutputs);

    PyObject *pMaxArgPos = PyObject_CallMethod
    (
        self->pModel, "inputs_max_argPos", NULL
    );
    if(!pMaxArgPos){ Modena_PyErr_Print(); }
    self->inputs_size = 1 + PyInt_AsSsize_t(pMaxArgPos);
    Py_DECREF(pMaxArgPos);

    self->inherited_inputs_size = 0;

    // Avoiding double indirection in modena_model_call
    // Use modena_function_new to construct, then copy function pointer
    self->mf = modena_function_new(self);
    self->function = self->mf->function;

    self->argPos_used = malloc
    (
        (self->inputs_size + self->inherited_inputs_size)*sizeof(bool)
    );

    modena_model_read_substituteModels(self);

    if(!pParameters)
    {
        pParameters = PyObject_GetAttrString(self->pModel, "parameters");
        if(!pParameters){ Modena_PyErr_Print(); }
    }
    else
    {
        Py_INCREF(pParameters);
    }

    PyObject *pSeq = PySequence_Fast(pParameters, "expected a sequence");
    self->parameters_size = PySequence_Size(pParameters);
    self->parameters = malloc(self->parameters_size*sizeof(double));
    size_t i;
    for(i = 0; i < self->parameters_size; i++)
    {
        self->parameters[i] = PyFloat_AsDouble(PyList_GET_ITEM(pSeq, i));
    }
    Py_DECREF(pSeq);
    Py_DECREF(pParameters);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    return 0;
}

static PyMemberDef modena_model_t_members[] = {
    {NULL}  /* Sentinel */
};

static PyObject *
modena_model_t_call(modena_model_t* self, PyObject *args, PyObject *kwds)
{
    PyObject *pIn_i=NULL, *pI=NULL, *pCheckBounds=NULL;
    bool checkBounds = true;

    static char *kwlist[] = {"inputs", "inherited_inputs", "checkBounds"};

    if
    (
        !PyArg_ParseTupleAndKeywords
        (
            args, kwds, "OO|O", kwlist, &pIn_i, &pI, &pCheckBounds
        )
    )
    {
        printf("Expected two arguments\n");
        return NULL;
    }

    if(pCheckBounds)
    {
        checkBounds = PyObject_IsTrue(pCheckBounds);
    }

    if(!PyList_Check(pIn_i))
    {
        printf("First argument is not a list\n");
        return NULL;
    }

    if(!PyList_Check(pI))
    {
        printf("First argument is not a list\n");
        return NULL;
    }

    modena_inputs_t *inputs = modena_inputs_new(self);

    PyObject *pSeq = PySequence_Fast(pIn_i, "expected a sequence");
    size_t len = PySequence_Size(pIn_i);
    size_t j;
    for(j = 0; j < len; j++)
    {
        modena_inherited_inputs_set
        (
            inputs, j, PyFloat_AsDouble(PyList_GET_ITEM(pSeq, j))
        );
    }
    Py_DECREF(pSeq);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    pSeq = PySequence_Fast(pI, "expected a sequence");
    len = PySequence_Size(pI);
    for(j = 0; j < len; j++)
    {
        modena_inputs_set
        (
            inputs, j, PyFloat_AsDouble(PyList_GET_ITEM(pSeq, j))
        );
    }
    Py_DECREF(pSeq);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    modena_outputs_t *outputs = modena_outputs_new(self);

    if(checkBounds)
    {
        modena_model_call(self, inputs, outputs);
    }
    else
    {
        modena_model_call_no_check(self, inputs, outputs);
    }

    PyObject* pOutputs = PyList_New(self->outputs_size);
    for(j = 0; j < self->outputs_size; j++)
    {
        PyList_SET_ITEM
        (
            pOutputs, j, PyFloat_FromDouble(modena_outputs_get(outputs, j))
        );
    }
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    modena_inputs_destroy(inputs);
    modena_outputs_destroy(outputs);

    return pOutputs;
}

static PyMethodDef modena_model_t_methods[] = {
    {"call", (PyCFunction) modena_model_t_call, METH_KEYWORDS,
     "Return the name, combining the first and last name"
    },
    {NULL}  /* Sentinel */
};

static PyTypeObject modena_model_tType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "modena.modena_model_t", /*tp_name*/
    sizeof(modena_model_t), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)modena_model_t_dealloc, /*tp_dealloc*/
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
    "modena_model_t objects", /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    modena_model_t_methods, /* tp_methods */
    modena_model_t_members, /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)modena_model_t_init, /* tp_init */
    0,                         /* tp_alloc */
    modena_model_t_new,  /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initlibmodena(void)
{
    PyObject* m;

    if(PyType_Ready(&modena_model_tType) < 0)
    {
        return;
    }

    m = Py_InitModule3
    (
        "libmodena",
        module_methods,
        "Example module that creates an extension type."
    );

    if(m == NULL)
    {
        return;
    }

    Py_INCREF(&modena_model_tType);
    PyModule_AddObject(m, "modena_model_t", (PyObject *) &modena_model_tType);
}

// modena_model_t instantiates the surrogate model
modena_model_t *modena_model_new
(
    const char *modelId
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
    PyObject *kw = Py_BuildValue("{s:s}", "modelId", modelId);

    PyObject *m = PyObject_Call((PyObject *) &modena_model_tType, args, kw);
    Py_DECREF(args);
    Py_DECREF(kw);
    if(!m){ Modena_PyErr_Print(); }

    return (modena_model_t *) m;
}

size_t modena_model_inputs_argPos(const modena_model_t *m, const char *name)
{
    PyObject *pRet = PyObject_CallMethod
    (
        m->pModel, "inputs_argPos", "(z)", name
    );
    if(!pRet){ Modena_PyErr_Print(); }
    size_t argPos = PyInt_AsSsize_t(pRet);
    Py_DECREF(pRet);

    if(m->argPos_used)
    {
        m->argPos_used[argPos] = true;
    }

    return argPos;
}

size_t modena_model_inherited_inputs_argPos
(
    const modena_model_t *m,
    const char *name
)
{
    fprintf(stderr, "Not implemented\n");
    exit(1);
}

size_t modena_model_outputs_argPos(const modena_model_t *m, const char *name)
{
    PyObject *pRet = PyObject_CallMethod
    (
        m->pModel, "outputs_argPos", "(z)", name
    );
    if(!pRet){ Modena_PyErr_Print(); }
    size_t ret = PyInt_AsSsize_t(pRet);
    Py_DECREF(pRet);

    return ret;
}

void modena_model_argPos_check(const modena_model_t *m)
{
    bool allUsed = true;
    size_t i = 0;
    size_t j = 0;

    for(j = 0; j < m->inputs_size; j++)
    {
        if(!m->argPos_used[i++])
        {
            //TODO: Replace by call into python
            //printf("argPos for %s not used\n", m->inputs_names[j]);
            allUsed = false;
        }
    }

    if(!allUsed)
    {
        fprintf(stderr, "Not all input arguments used - Exiting\n");
        exit(1);
    }

    for(j = 0; j < m->inherited_inputs_size; j++)
    {
        if(!m->argPos_used[i++])
        {
            allUsed = false;
        }
    }

    if(!allUsed)
    {
        fprintf(stderr, "Not all inherited input arguments used\n");
        exit(1);
    }
}

size_t modena_model_inputs_size(const modena_model_t *m)
{
    return m->inputs_size;
}

size_t modena_model_inherited_inputs_size(const modena_model_t *m)
{
    return m->inherited_inputs_size;
}

size_t modena_model_outputs_size(const modena_model_t *m)
{
    return m->outputs_size;
}

void modena_model_inputs_siunits
(
    const modena_model_t *m,
    const size_t i,
    modena_siunits_t *units
)
{
    fprintf(stderr, "Not implemented\n");
    exit(1);
}

void modena_model_inherited_inputs_siunits
(
    const modena_model_t *m,
    const size_t i,
    modena_siunits_t *units
)
{
    fprintf(stderr, "Not implemented\n");
    exit(1);
}

void modena_model_outputs_siunits
(
    const modena_model_t *m,
    const size_t i,
    modena_siunits_t *units
)
{
    fprintf(stderr, "Not implemented\n");
    exit(1);
}

void modena_substitute_model_call
(
    const modena_substitute_model_t *sm,
    const modena_model_t *parent,
    modena_inputs_t *inputs
)
{
    size_t j;
    for(j = 0; j < sm->map_inputs_size; j++)
    {
        //printf("i%zu <- ip%zu\n", sm->map_inputs[2*j+1], sm->map_inputs[2*j]);
        sm->inputs->inputs[sm->map_inputs[2*j+1]] =
            inputs->inputs[sm->map_inputs[2*j]];
    }

    modena_model_call(sm->model, sm->inputs, sm->outputs);

    for(j = 0; j < sm->map_outputs_size; j++)
    {
        //printf("ip%zu <- o%zu\n", sm->map_outputs[2*j+1], sm->map_outputs[2*j]);
        inputs->inputs[sm->map_outputs[2*j+1]] =
            sm->outputs->outputs[sm->map_outputs[2*j]];
    }
}

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
    modena_model_t *m,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
)
{
    size_t j;
    for(j = 0; j < m->substituteModels_size; j++)
    {
        modena_substitute_model_call(&m->substituteModels[j], m, inputs);
    }

    for(j = 0; j < m->inputs_minMax_size; j++)
    {
        /*
        printf
        (
            "j = %zu %g < %g || %g > %g\n",
            j,
            inputs->inputs[j],
            m->inputs_min[j],
            inputs->inputs[j],
            m->inputs_max[j]
        );
        */

        if
        (
            inputs->inputs[j] < m->inputs_min[j]
         || inputs->inputs[j] > m->inputs_max[j]
        )
        {
            PyObject* pOutside = PyList_New(m->inputs_minMax_size);

            for(j = 0; j < m->inputs_minMax_size; j++)
            {
                PyList_SET_ITEM
                (
                    pOutside, j, PyFloat_FromDouble(inputs->inputs[j])
                );
            }

            PyObject *pRet = PyObject_CallMethod
            (
                m->pModel, "outOfBounds", "(O)", pOutside
            );
            Py_DECREF(pOutside);
            if(!pRet){ Modena_PyErr_Print(); }
            int ret = PyInt_AsLong(pRet);
            Py_DECREF(pRet);

            return ret;
        }
    }

    m->function
    (
        m->parameters,
        inputs->inherited_inputs,
        inputs->inputs,
        outputs->outputs
    );

    return 0;
}

void modena_model_call_no_check
(
    modena_model_t *m,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
)
{
    size_t j;
    for(j = 0; j < m->substituteModels_size; j++)
    {
        modena_substitute_model_call(&m->substituteModels[j], m, inputs);
    }

    m->function
    (
        m->parameters,
        inputs->inherited_inputs,
        inputs->inputs,
        outputs->outputs
    );
}

void modena_model_destroy(modena_model_t *m)
{
    size_t j;
    for(j = 0; j < m->substituteModels_size; j++)
    {
        Py_DECREF(m->substituteModels[j].model);
        modena_inputs_destroy(m->substituteModels[j].inputs);
        modena_outputs_destroy(m->substituteModels[j].outputs);
        free(m->substituteModels[j].map_inputs);
        free(m->substituteModels[j].map_outputs);
    }
    free(m->substituteModels);

    free(m->parameters);
    free(m->inputs_min);
    free(m->inputs_max);

    free(m->argPos_used);

    modena_function_destroy(m->mf);

    Py_DECREF(m->pModel);

    m->ob_type->tp_free((PyObject*)m);
}

// modena_model_t instantiates a surrogate function
modena_function_t *modena_function_new
(
    const modena_model_t *m
)
{
    modena_function_t *f = malloc(sizeof(modena_function_t));

    if(lt_dlinit())
    {
        fprintf(stderr, "lt_dlinit: %s\n", lt_dlerror());
        exit(1);
    }

    PyObject *pSurrogateFunction =
        PyObject_GetAttrString(m->pModel, "surrogateFunction");
    if(!pSurrogateFunction){ Modena_PyErr_Print(); }

    PyObject *pFunctionName =
        PyObject_GetAttrString(pSurrogateFunction, "functionName");
    if(!pFunctionName){ Modena_PyErr_Print(); }

    PyObject *pLibraryName =
        PyObject_GetAttrString(pSurrogateFunction, "libraryName");
    if(!pLibraryName){ Modena_PyErr_Print(); }

    Py_DECREF(pSurrogateFunction);

    f->handle = lt_dlopen(PyString_AsString(pLibraryName));
    if(!f->handle)
    {
        fprintf
        (
           stderr,
           "lt_dlopen: Could not open library %s\nlt_dlopen: %s\n",
           PyString_AsString(pLibraryName),
           lt_dlerror()
        );
        exit(1);
    }

    f->function = lt_dlsym(f->handle, PyString_AsString(pFunctionName));
    if(!f->function)
    {
        fprintf
        (
            stderr,
            "lt_dlsym: Could not find function %s in library %s\n"
            "lt_dlsym: %s",
            PyString_AsString(pFunctionName),
            PyString_AsString(pLibraryName),
            lt_dlerror()
        );
        lt_dlclose(f->handle);
        exit(1);
    }

    Py_DECREF(pFunctionName);
    Py_DECREF(pLibraryName);

    return f;
}

void modena_function_destroy(modena_function_t *f)
{
    lt_dlclose(f->handle);
    lt_dlexit();
    free(f);
}

