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

#include "model.h"
#include "structmember.h"
#include "global.h"

PyObject *modena_SurrogateModel = NULL;

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

bool modena_model_read_substituteModels(modena_model_t *self)
{
    //Modena_Info_Print("In %s", __func__);

    PyObject *pSubstituteModels = PyObject_GetAttrString
    (
        self->pModel, "substituteModels"
    );
    if(!pSubstituteModels){ Modena_PyErr_Print(); }

    PyObject *pSeq = PySequence_Fast
    (
        pSubstituteModels, "expected a sequence"
    );
    self->substituteModels_size = PySequence_Size(pSubstituteModels);
    self->substituteModels =
        malloc(self->substituteModels_size*sizeof(modena_substitute_model_t));
    size_t i;
    for(i = 0; i < self->substituteModels_size; i++)
    {
        PyObject *args = PyTuple_New(0);
        PyObject *kw = Py_BuildValue
        (
            "{s:O}", "model", PyList_GET_ITEM(pSeq, i)
        );

        self->substituteModels[i].model = (modena_model_t *) PyObject_Call
        (
            (PyObject *) &modena_model_tType,
            args,
            kw
        );
        Py_DECREF(args);
        Py_DECREF(kw);

        if(!self->substituteModels[i].model)
        {
            if
            (
                PyErr_ExceptionMatches(modena_DoesNotExist)
             || PyErr_ExceptionMatches(modena_ParametersNotValid)
            )
            {
                PyObject *pModelId =
                    PyObject_GetAttrString(PyList_GET_ITEM(pSeq, i), "_id");
                if(!pModelId){ Modena_PyErr_Print(); }
                const char* modelId = PyString_AsString(pModelId);
                Py_DECREF(pModelId);

                PyObject *pRet = NULL;
                if
                (
                    PyErr_ExceptionMatches(modena_DoesNotExist)
                )
                {
                    fprintf
                    (
                        stderr,
                        "Loading model %s failed - Attempting automatic initialisation\n",
                        modelId
                    );

                    pRet = PyObject_CallMethod
                    (
                        modena_SurrogateModel,
                        "exceptionLoad",
                        "(z)",
                        modelId
                    );
                }
                else
                {
                    fprintf
                    (
                        stderr,
                        "Parameters of model %s are invalid - Trying to initialise\n",
                        modelId
                    );

                    pRet = PyObject_CallMethod
                    (
                        modena_SurrogateModel,
                        "exceptionParametersNotValid",
                        "(z)",
                        modelId
                    );
                }

                if(!pRet){ Modena_PyErr_Print(); }
                int ret = PyInt_AsLong(pRet);
                Py_DECREF(pRet);

                modena_error_code = ret;

                Py_DECREF(pSeq);
                Py_DECREF(pSubstituteModels);

                return false;
            }
            else
            {
                Modena_PyErr_Print();
                return false;
            }
        }

        self->substituteModels[i].inputs = modena_inputs_new
        (
            self->substituteModels[i].model
        );

        self->substituteModels[i].outputs = modena_outputs_new
        (
            self->substituteModels[i].model
        );

        modena_substitute_model_calculate_maps
        (
            &self->substituteModels[i],
            self
        );
    }

    Py_DECREF(pSeq);
    Py_DECREF(pSubstituteModels);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    return true;
}

void modena_model_get_minMax
(
    modena_model_t *self
)
{
    PyObject *pObj = PyObject_CallMethod(self->pModel, "minMax", NULL);
    if(!pObj){ Modena_PyErr_Print(); }

    PyObject *pMin = PyTuple_GET_ITEM(pObj, 0); // Borrowed ref
    PyObject *pSeq = PySequence_Fast(pMin, "expected a sequence");
    self->inputs_size = PySequence_Size(pMin);
    self->inputs_min = malloc(self->inputs_size*sizeof(double));
    size_t i;
    for(i = 0; i < self->inputs_size; i++)
    {
        self->inputs_min[i] = PyFloat_AsDouble(PyList_GET_ITEM(pSeq, i));
    }
    Py_DECREF(pSeq);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    PyObject *pMax = PyTuple_GET_ITEM(pObj, 1); // Borrowed ref
    pSeq = PySequence_Fast(pMax, "expected a sequence");
    self->inputs_max = malloc(self->inputs_size*sizeof(double));
    for(i = 0; i < self->inputs_size; i++)
    {
        self->inputs_max[i] = PyFloat_AsDouble(PyList_GET_ITEM(pSeq, i));
    }
    Py_DECREF(pSeq);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    Py_DECREF(pObj);
}

modena_model_t *modena_model_new
(
    const char *modelId
)
{
    PyObject *args = PyTuple_New(0);
    PyObject *kw = Py_BuildValue("{s:s}", "modelId", modelId);

    PyObject *pNewObj = PyObject_Call
    (
        (PyObject *) &modena_model_tType,
        args,
        kw
    );

    Py_DECREF(args);
    Py_DECREF(kw);
    if(!pNewObj)
    {
        if
        (
            PyErr_ExceptionMatches(modena_DoesNotExist)
         || PyErr_ExceptionMatches(modena_ParametersNotValid)
        )
        {
            PyErr_Clear();

            PyObject *pRet = NULL;
            if
            (
                PyErr_ExceptionMatches(modena_DoesNotExist)
            )
            {
                fprintf
                (
                    stderr,
                    "Loading model %s failed - "
                    "Attempting automatic initialisation\n",
                    modelId
                );

                pRet = PyObject_CallMethod
                (
                    modena_SurrogateModel,
                    "exceptionLoad",
                    "(z)",
                    modelId
                );
            }
            else
            {
                fprintf
                (
                    stderr,
                    "Parameters of model %s are invalid - "
                    "Trying to initialise\n",
                    modelId
                );

                pRet = PyObject_CallMethod
                (
                    modena_SurrogateModel,
                    "exceptionParametersNotValid",
                    "(z)",
                    modelId
                );
            }

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

    return (modena_model_t *) pNewObj;
}

size_t modena_model_inputs_argPos(const modena_model_t *self, const char *name)
{
    PyObject *pRet = PyObject_CallMethod
    (
        self->pModel,
        "inputs_argPos",
        "(z)",
        name
    );
    if(!pRet){ Modena_PyErr_Print(); }
    size_t argPos = PyInt_AsSsize_t(pRet);
    Py_DECREF(pRet);

    if(self->argPos_used)
    {
        //Modena_Info_Print
        //(
        //    "Mark argPos %zu as used from inputs_argPos\n",
        //    argPos
        //);
        self->argPos_used[argPos] = true;
    }

    return argPos;
}

size_t modena_model_outputs_argPos(const modena_model_t *self, const char *name)
{
    PyObject *pRet = PyObject_CallMethod
    (
        self->pModel,
        "outputs_argPos",
        "(z)",
        name
    );
    if(!pRet){ Modena_PyErr_Print(); }
    size_t ret = PyInt_AsSsize_t(pRet);
    Py_DECREF(pRet);

    return ret;
}

void modena_model_argPos_check(const modena_model_t *self)
{
    bool allUsed = true;
    size_t j = 0;

    for(j = 0; j < self->inputs_size; j++)
    {
        if(!self->argPos_used[j])
        {
            //TODO: Replace by call into python
            //Modena_Info_Print("argPos for %s not used", self->inputs_names[j]);
            fprintf(stderr, "argPos %zu not used", j);
            allUsed = false;
            break;
        }
    }

    if(!allUsed)
    {
        fprintf(stderr, "Not all input arguments used - Exiting\n");
        exit(1);
    }
}

size_t modena_model_inputs_size(const modena_model_t *self)
{
    return self->inputs_size;
}

size_t modena_model_outputs_size(const modena_model_t *self)
{
    return self->outputs_size;
}

int modena_substitute_model_call
(
    const modena_substitute_model_t *sm,
    const modena_model_t *parent,
    modena_inputs_t *inputs
)
{
    size_t j;
    for(j = 0; j < sm->map_inputs_size; j++)
    {
        /*
        printf
        (
            "i%zu <- ip%zu (%g)\n",
            sm->map_inputs[2*j+1],
            sm->map_inputs[2*j],
            inputs->inputs[sm->map_inputs[2*j]]
        );
        */
        sm->inputs->inputs[sm->map_inputs[2*j+1]] =
            inputs->inputs[sm->map_inputs[2*j]];
    }

    int ret = modena_model_call(sm->model, sm->inputs, sm->outputs);
    if(ret){ return ret; }

    for(j = 0; j < sm->map_outputs_size; j++)
    {
        /*
        printf
        (
            "ip%zu <- o%zu (%g)\n",
            sm->map_outputs[2*j+1],
            sm->map_outputs[2*j],
            sm->outputs->outputs[sm->map_outputs[2*j]]
        );
        */
        inputs->inputs[sm->map_outputs[2*j+1]] =
            sm->outputs->outputs[sm->map_outputs[2*j]];
    }

    return 0;
}

int write_outside_point
(
    modena_model_t *self,
    modena_inputs_t *inputs
)
{
    PyObject* pOutside = PyList_New(self->inputs_size);

    size_t j;
    for(j = 0; j < self->inputs_size; j++)
    {
        PyList_SET_ITEM
        (
            pOutside, j, PyFloat_FromDouble(inputs->inputs[j])
        );
    }

    PyObject *pRet = PyObject_CallMethod
    (
       self->pModel,
       "exceptionOutOfBounds",
       "(O)",
       pOutside
    );
    Py_DECREF(pOutside);
    if(!pRet){ Modena_PyErr_Print(); }
    int ret = PyInt_AsLong(pRet);
    Py_DECREF(pRet);

    modena_error_code = ret;

    return ret;
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
    modena_model_t *self,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
)
{
    if
    (
          self->parameters_size == 0
       && self->parameters_size != self->mf->parameters_size
    )
    {
        return write_outside_point(self, inputs);
    }

    size_t j;
    for(j = 0; j < self->substituteModels_size; j++)
    {
        int ret = modena_substitute_model_call
        (
            &self->substituteModels[j],
            self,
            inputs
        );
        if(ret){ return ret; }
    }

    for(j = 0; j < self->inputs_size; j++)
    {
        /*
        printf
        (
            "j = %zu %g < %g || %g > %g\n",
            j,
            inputs->inputs[j],
            self->inputs_min[j],
            inputs->inputs[j],
            self->inputs_max[j]
        );
        */

        if
        (
            inputs->inputs[j] < self->inputs_min[j]
         || inputs->inputs[j] > self->inputs_max[j]
        )
        {
            return write_outside_point(self, inputs);
        }
    }

    self->function
    (
        self,
        inputs->inputs,
        outputs->outputs
    );

    return 0;
}

void modena_model_call_no_check
(
    modena_model_t *self,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
)
{
    //Modena_Info_Print("In %s", __func__);

    if
    (
          self->parameters_size == 0
       && self->parameters_size != self->mf->parameters_size
    )
    {
        write_outside_point(self, inputs);
    }

    size_t j;
    for(j = 0; j < self->substituteModels_size; j++)
    {
        modena_substitute_model_call
        (
            &self->substituteModels[j],
            self,
            inputs
        );
    }

    for(j = 0; j < self->inputs_size; j++)
    {
        /*
        printf
        (
            "j = %zu %g\n",
            j,
            inputs->inputs[j]
        );
        */
    }

    self->function
    (
        self,
        inputs->inputs,
        outputs->outputs
    );
}

/* Destructor, frees the memory block occupied by a model.
 */
void modena_model_destroy(modena_model_t *self)
{
    size_t j;
    for(j = 0; j < self->substituteModels_size; j++)
    {
        Py_XDECREF(self->substituteModels[j].model);
        modena_inputs_destroy(self->substituteModels[j].inputs);
        modena_outputs_destroy(self->substituteModels[j].outputs);
        free(self->substituteModels[j].map_inputs);
        free(self->substituteModels[j].map_outputs);
    }
    free(self->substituteModels);

    free(self->parameters);
    free(self->inputs_min);
    free(self->inputs_max);

    free(self->argPos_used);

    if(self->mf)
    {
        modena_function_destroy(self->mf);
    }

    Py_XDECREF(self->pModel);

    self->ob_type->tp_free((PyObject*)self);
}

/* C-Python: Destructor, exposed as __del__ in Python
 */
static void modena_model_t_dealloc(modena_model_t* self)
{
    modena_model_destroy(self);
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
static PyMemberDef modena_model_t_members[] = {
    {"outputs_size", T_PYSSIZET,
       offsetof(modena_model_t, outputs_size), READONLY , "number of putputs"},
    {"inputs_size", T_PYSSIZET,
      offsetof(modena_model_t, inputs_size), READONLY , "number of inputs"},
    {"parameters_size", T_PYSSIZET,
      offsetof(modena_model_t, parameters_size), READONLY , "number of parameters"},
    {NULL}  /* Sentinel */
};

/* C-Python: Method exposed in Python as __call__
 *
 * TODO: The method is also exposed as "call", but this should be deprecated
 */
static PyObject *modena_model_t_call
(
    modena_model_t* self,
    PyObject *args,
    PyObject *kwds
)
{
    //Modena_Info_Print("In %s", __func__);

    PyObject *pI=NULL, *pCheckBounds=NULL;
    bool checkBounds = true;

    static char *kwlist[] = { "inputs", "checkBounds", NULL };

    if
    (
        !PyArg_ParseTupleAndKeywords
        (
            args,
            kwds,
            "O|O",
            kwlist,
            &pI,
            &pCheckBounds
        )
    )
    {
        Modena_PyErr_Print();
    }

    if(pCheckBounds)
    {
        checkBounds = PyObject_IsTrue(pCheckBounds);
    }

    if(!PyList_Check(pI))
    {
        printf("First argument is not a list\n");
        return NULL;
    }

    PyObject *pSeq = PySequence_Fast(pI, "expected a sequence");
    size_t len = PySequence_Size(pI);

    if(len != self->inputs_size)
    {
        Py_DECREF(pSeq);
        printf("input array has incorrect size %zu %zu\n", len, self->inputs_size);
        return NULL;
    }

    modena_inputs_t *inputs = modena_inputs_new(self);

    size_t j;
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
        if(modena_model_call(self, inputs, outputs))
        {
            modena_inputs_destroy(inputs);
            modena_outputs_destroy(outputs);

            PyErr_SetString
            (
                modena_OutOfBounds,
                "Surrogate model is used out-of-bounds"
            );

            return NULL;
        }
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
static PyMethodDef modena_model_t_methods[] = {
    {"call", (PyCFunction) modena_model_t_call, METH_KEYWORDS,
        "Call surrogate model and return outputs"
    },
    {NULL}  /* Sentinel */
};

/*
 */
PyObject*
modena_model_t_get_parameters(modena_model_t *self, void *closure)
{
    PyObject* pParams = PyList_New(self->parameters_size);
    size_t i;
    for(i = 0; i < self->parameters_size; i++)
    {
        PyList_SET_ITEM(pParams, i, PyFloat_FromDouble(self->parameters[i]) );
    }
    return pParams;
}

/*
 */
static int
modena_model_t_set_parameters(modena_model_t *self, PyObject *value, void *closure)
{
   // TODO: Error checks for the following cases:
   //       1. len(value) == self->parameters_size
   //       2. type(value) == list or tuple
   //       3. value != NULL

    /*if (value == NULL)
    {
          PyErr_SetString(PyExc_TypeError, "Cannot delete parameter values");
          return -1;
    }
    if (! PyString_Check(value)) {
          PyErr_SetString(PyErr_TypeError, "First attribute must be a string");
          return -1;
    }*/

    size_t i;
    for(i = 0; i < self->parameters_size; i++)
    {
         self->parameters[i] = PyFloat_AsDouble(PyList_GetItem(value, i));
    }

    // PyErr_SetString(PyExc_TypeError, "Attribute is read-only!");
    return 0;
}

/* C-Python
 */
static PyGetSetDef modena_model_t_getset[] = {
    {"parameters",
      (getter)modena_model_t_get_parameters,
      (setter)modena_model_t_set_parameters,
     "parameters",
      NULL},
    {NULL} /* Sentinel */
};

/* C-Python: Initialiser, exposed in Python as the method: __new__
 */
static int modena_model_t_init
(
    modena_model_t *self,
    PyObject *args,
    PyObject *kwds
)
{
    //Modena_Info_Print("In %s", __func__);

    PyObject *pParameters=NULL, *pModel=NULL;
    char *modelId=NULL;
    size_t i, j;

    static char *kwlist[] = {"model", "modelId", "parameters", NULL};

    if
    (
        !PyArg_ParseTupleAndKeywords
        (
            args,
            kwds,
            "|OsO",
            kwlist,
            &pModel,
            &modelId,
            &pParameters
        )
    )
    {
        Modena_PyErr_Print();
    }

    if(!pModel)
    {
        self->pModel = PyObject_CallMethod
        (
            modena_SurrogateModel,
            "load",
            "(z)",
            modelId
        );

        if(!self->pModel)
        {
            PyErr_SetString
            (
                modena_DoesNotExist,
                "Surrogate model does not exist"
            );

            return -1;
        }
    }
    else
    {
        Py_INCREF(pModel);
        self->pModel = pModel;
    }

    modena_model_get_minMax(self);

    //PyObject_Print(self->pModel, stdout, 0);
    //printf("\n");

    PyObject *pOutputs = PyObject_GetAttrString(self->pModel, "outputs");
    if(!pOutputs){ Modena_PyErr_Print(); }
    self->outputs_size = PyDict_Size(pOutputs);
    Py_DECREF(pOutputs);

    if(!modena_model_read_substituteModels(self))
    {
        return -1;
    }

    // Avoiding double indirection in modena_model_call
    // Use modena_function_new to construct, then copy function pointer
    self->mf = modena_function_new_from_model(self);
    self->function = self->mf->function;

    self->argPos_used = malloc(self->inputs_size*sizeof(bool));

    for(j = 0; j < self->inputs_size; j++)
    {
        self->argPos_used[j] = false;
    }

    for(j = 0; j < self->substituteModels_size; j++)
    {
        modena_substitute_model_t *sm = &self->substituteModels[j];
        for(i = 0; i < sm->map_outputs_size; i++)
        {
            //printf("Mark argPos %zu as used\n", sm->map_outputs[2*i+1]);
            self->argPos_used[sm->map_outputs[2*i+1]] = true;
        }
    }

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

    if
    (
          self->parameters_size == 0
       && self->parameters_size != self->mf->parameters_size
    )
    {
        PyObject *args = PyTuple_New(2);
        PyObject* str = PyString_FromString
        (
            "Surrogate model does not have valid parameters"
        );
        PyTuple_SET_ITEM(args, 0, str);
        PyTuple_SET_ITEM(args, 1, self->pModel);

        PyErr_SetObject
        (
            modena_ParametersNotValid,
            args
        );

        Py_DECREF(pSeq);
        Py_DECREF(pParameters);
        return -1;
    }

    self->parameters = malloc(self->parameters_size*sizeof(double));
    for(i = 0; i < self->parameters_size; i++)
    {
        self->parameters[i] = PyFloat_AsDouble(PyList_GET_ITEM(pSeq, i));
    }
    Py_DECREF(pSeq);
    Py_DECREF(pParameters);
    if(PyErr_Occurred()){ Modena_PyErr_Print(); }

    return 0;
}

/* C-Python: Constructor, exposed in Python as the method: __new__
 */
static PyObject * modena_model_t_new
(
    PyTypeObject *type,
    PyObject *args,
    PyObject *kwds
)
{
    modena_model_t *self;

    self = (modena_model_t *)type->tp_alloc(type, 0);
    if(self)
    {
        // Set everything to zero
        self->pModel = NULL;
        self->outputs_size = 0;
        self->inputs_size = 0;
        self->inputs_min = NULL;
        self->inputs_max = NULL;
        self->argPos_used = NULL;
        self->parameters_size = 0;
        self->parameters = NULL;
        self->mf = NULL;
        self->function = NULL;
        self->substituteModels_size = 0;
        self->substituteModels = NULL;
    }

    return (PyObject *)self;
}

/* C-Python: The C structure used to describe the modena_model type.
 */
PyTypeObject modena_model_tType = {
    PyObject_HEAD_INIT(NULL)
    0,                                                              /*ob_size*/
    "modena.modena_model_t",                                        /*tp_name*/
    sizeof(modena_model_t),                                    /*tp_basicsize*/
    0,                                                          /*tp_itemsize*/
    (destructor)modena_model_t_dealloc,                          /*tp_dealloc*/
    0,                                                             /*tp_print*/
    0,                                                           /*tp_getattr*/
    0,                                                           /*tp_setattr*/
    0,                                                           /*tp_compare*/
    0,                                                              /*tp_repr*/
    0,                                                         /*tp_as_number*/
    0,                                                       /*tp_as_sequence*/
    0,                                                        /*tp_as_mapping*/
    0,                                                             /*tp_hash */
    (ternaryfunc)modena_model_t_call,                               /*tp_call*/
    0,                                                               /*tp_str*/
    0,                                                          /*tp_getattro*/
    0,                                                          /*tp_setattro*/
    0,                                                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                      /*tp_flags*/
    "modena_model_t objects",                                      /* tp_doc */
    0,                                                        /* tp_traverse */
    0,                                                           /* tp_clear */
    0,                                                     /* tp_richcompare */
    0,                                                  /* tp_weaklistoffset */
    0,                                                            /* tp_iter */
    0,                                                        /* tp_iternext */
    modena_model_t_methods,                                    /* tp_methods */
    modena_model_t_members,                                    /* tp_members */
    modena_model_t_getset,                                      /* tp_getset */
    0,                                                            /* tp_base */
    0,                                                            /* tp_dict */
    0,                                                       /* tp_descr_get */
    0,                                                       /* tp_descr_set */
    0,                                                      /* tp_dictoffset */
    (initproc)modena_model_t_init,                                /* tp_init */
    0,                                                           /* tp_alloc */
    modena_model_t_new,                                            /* tp_new */
};

