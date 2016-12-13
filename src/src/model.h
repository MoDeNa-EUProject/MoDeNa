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

MoDeNa low-level interface library

@author Henrik Rusche
@copyright  2014-2016, MoDeNa Project. GNU Public License.
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
 * @brief stores a model and mapping for substitution
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
 * @brief stores a surrogate model
 *
*/
typedef struct modena_model_t
{
    PyObject_HEAD

    PyObject *pModel;                        /**< Reference to python object.*/

    size_t outputs_size;                        /**< Length of output vector.*/

    size_t inputs_size;                         /**< Length of input vector. */

    double *inputs_min;            /**< Lower bounds of the input arguments. */

    double *inputs_max;            /**< Upper bounds of the input arguments. */

    bool *argPos_used;         /**< Indicators used to track argument usage. */

    size_t parameters_size;                 /**< Length of parameter vector. */

    double *parameters;                /**< Surrogate model parameter values */

    struct modena_function_t *mf; /**< Surrogate function `modena_function_t`*/

    void (*function)
    (
        const struct modena_model_t* model,
        const double* i,
        double *o
    );               /**< Function pointer to pre-compiled surrogate function*/

    size_t substituteModels_size;    /**< Length of substitute models vector */

    modena_substitute_model_t *substituteModels; /**< Substitute models `modena_substitute_model_t`*/

} modena_model_t;

/**
 * @brief Function fetching a surrogate model from MongoDB.
 *
 * ### About writing adaptors for surrogate models
 *
 * A adaptor is a code fragment which makes an application able to use the
 * surrogate models (SM) that are stored in the MongoDB database.
 * Writing a adaptor for a SM requires implementation of
 * code fragments that corresponds to the life-span of a surrogate model, which
 * consists of three phases:
 *   * Initialisation
 *     1. Fetch the model from database.
 *     2. Allocate memory for input and output vectors.
 *     3. Query the SM for the position of each individual input and output.
 *   * Execution
 *     1. Set the input vector.
 *     2. Evaluate the surrogate model.
 *     3. Check the MoDeNa framework for errors.
 *     4. Fetch outputs.
 *   * Termination
 *     1. Deallocate memory.
 *
 * ### About
 *
 * The function `modena_model_new` is used in the initialisation phase of the
 * adaptor, and its purpose is to fetch a surrogate model from the database.
 * The input to the function is the name, technically the database "_id", of
 * the surrogate model.
 *
 * When the surrogate model has been fetched from the database the
 * initialisation continues with allocating memory for the input and output
 * vectors. However, this procedure is only performed one time for every
 * surrogate model.
 *
 * ### Usage
 *
 * The function is only called one time for every surrogate model that the user
 * want to employ in a application. It is implemented as a pointer to
 * `modena_model_t` as follows:
 *
 * * C:
 * ~~~~{.c}
 * modena_model_t *model = modena_model_new("MY_MODEL");
 * ~~~~
 * * Fortran:
 * ~~~~{.f90}
 * type(c_ptr) :: model = c_null_ptr
 * model = modena_model_new (c_char_"MY_MODEL"//c_null_char);
 * ~~~~
 * * Python:
 * ~~~~{.py}
 * model = SurrogateModel.load("MY_MODEL")
 * ~~~~
 *
 * #### Important
 *
 * 1. Make sure that the name of the surrogate model is spelled correctly, i.e.
 *    that it corresponds to the "_id" field in the definition of the SM.
 *    ~~~~{.py}
 *        m = BackwardMappingModel(
 *               _id= "MY_MODEL",
 *               surrogateFunction= f,
 *               exactTask= FlowRateExactSim(),
 *               substituteModels= [ ],
 *               initialisationStrategy= Strategy.InitialPoints(),
 *               outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(),
 *               parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(),
 *            )
 *    ~~~~
 *
 * 2. Ensure that the input and output variables are spelled correctly,
 *    according to the surrogate function corresponding to the surrogate model.
 *    ~~~~{.py}
 *        f = CFunction(
 *          Ccode= ''' C-code Omitted ''',
 *          # Global bounds for the function
 *          inputs={
 *              'T': { 'min': -298.15, 'max': 5000.0 },
 *              'P': { 'min':       0, 'max':  100.0 },
 *          },
 *          outputs={
 *               'N': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
 *          },
 *          parameters={
 *              'param1': { 'min': 0.0, 'max': 10.0, 'argPos': 1 },
 *          },
 *        )
 *    ~~~~
 *
 * 3. Check 1 and 2.
 *
 * ~~~~{.c}
 * modena_model_t *model = modena_model_new("MY_MODEL");    // Fetch "MY_MODEL"
 *
 * modena_inputs_t *inputs = modena_inputs_new(model);        // Allocate input
 * modena_outputs_t *outputs = modena_outputs_new(model);    // Allocate output
 *
 * size_t T_pos = modena_model_inputs_argPos(model, "T"); // Input position "T"
 * size_t P_pos = modena_model_inputs_argPos(model, "P"); // Input position "P"
 * size_t N_pos = modena_model_outputs_argPos(model,"N");// Output position "N"
 *
 * modena_model_argPos_check(model);   // Check all positions have been queried
 * ~~~~
 *
 * The name of the model, here "MY_MODEL", must correspond to the "_id"
 * field in the definition of the surrogate model, which is located in a Python
 * module.
 *
 * #### Common issues:
 *
 * - A common error is to start a simulation without the surrogate model being
 *   located in the database. Check this by executing the line below in a 
 *   terminal (replacing "MY_MODEL" with the name of your surrogate model).
 *
 *   ~~~~{.sh}
 *   mongo --eval 'db.surrogate_model.find({"_id":"MY_MODEL"}).forEach(printjson)'
 *   ~~~~
 *
 * ---
 *
 * @param modelId (char) database '_id' if the desired surrogate model.
 * @return modena_model_t pointer to a surrogate model.
*/
modena_model_t *modena_model_new
(
    const char *modelId
);

/**
 * @brief Function determining position of an argument in the input vector.
 *
 * The function is used to determine the position of an argument @p name in the
 * input vector.
 *
 * @param self pointer to surrogate model created by `modena_model_new`.
 * @param name (char) name of the argument whose position is sought.
 * @return size_t integer representing the position of the argument.
*/
size_t modena_model_inputs_argPos
(
    const modena_model_t *self,
    const char *name
);

/**
 *  @brief Function checking that the user has queried all input positions.
 *  @param struct modena_model_t pointer to a surrogate model created by modena_model_new.
 *  @return void
*/
void modena_model_argPos_check(const modena_model_t *self);

/**
 *  @brief Function determining position of a result in the output vector.
 *
 *The function is used to determine the position of a result @p name in the
 *input vector.
 *
 *  @param self pointer to surrogate model created by modena_model_new.
 *  @param name (char) name of the result whose position is sought.
 *  @return size_t integer representing the position of the result.
*/
size_t modena_model_outputs_argPos
(
    const modena_model_t *self,
    const char *name
);

/**
 *  @brief Function returning the size of the input vector.
 *  @param modena_model_t pointer to a surrogate model created by modena_model_new.
 *  @return size_t integer length of the input array.
*/
size_t modena_model_inputs_size(const modena_model_t *self);

/**
 *  @brief Function returning the size of the output vector.
 *  @param modena_model_t surrogate model created by modena_model_new.
 *  @return size_t integer length of the output array.
*/
size_t modena_model_outputs_size(const modena_model_t *self);

void modena_model_inputs_siunits
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
 *  @brief Function calling the surrogate model and checking for errors.
 *  @param modena_model_t model pointer to a surrogate model.
 *  @param modena_inputs_t inputs pointer to the input vector
 *  @param modena_outputs_t outputs pointer to the output vector
 *  @return 201 requesting exit for new DOE without Restart
 *  @return 200 requesting exit for new DOE with Restart
 *  @return 100 updated model parameters, requesting to continue this run
 *  @return   1 failure
 *  @return   0  okay
 * If exit is requested, do what's necessary and exit with the same error code!
*/
int modena_model_call
(
    modena_model_t *model,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
);

/**
 *  @brief Function calling the surrogate model w/o checking for errors.
 *  @param model modena_model_t pointer to a surrogate model.
 *  @param inputs modena_inputs_t pointer to the input vector
 *  @param outputs modena_outputs_t pointer to the output vector
 *  @return void
*/
void modena_model_call_no_check
(
    modena_model_t *model,
    modena_inputs_t *inputs,
    modena_outputs_t *outputs
);

/**
 *  @brief Function deallocating the memory allocated for the surrogate model.
 *  @param model modena_model_t pointer to a surrogate model.
 *  @return void
*/
void modena_model_destroy(modena_model_t *model);

/** @} */ // end of C_interface_library

__END_DECLS

#endif /* __MODEL_H__ */

