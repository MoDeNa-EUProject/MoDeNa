/** @file modenaCalls.h
	@brief Allocate memory for the surrogate models
*/
#ifndef MODENADATA_H
#define MODENADATA_H

modena_model_t *bblgr1;
modena_model_t *bblgr2;

modena_model_t *kinetics;

modena_model_t *density_reaction_mixturemodel;

modena_inputs_t *inputs_bblgr1;
modena_outputs_t *outputs_bblgr1;

modena_inputs_t *inputs_bblgr2;
modena_outputs_t *outputs_bblgr2;

modena_inputs_t *inputs_kinetics;
modena_outputs_t *outputs_kinetics;

modena_inputs_t *inputs_den;
modena_outputs_t *outputs_den;

#endif
