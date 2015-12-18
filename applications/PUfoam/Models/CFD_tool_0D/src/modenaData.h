/** @file modenaCalls.h
	@brief instantiates the surrogate models, allocates memory and fetches arg positions
*/
#ifndef MODENADATA_H
#define MODENADATA_H
/*
instantiate the surrogate models:
    - bubbleGrowth1,
    - bubbleGrowth2,
    - density_reaction_mixture,
    - rheology,
    - simpleKinetics
*/

modena_model_t *bblgr1;
modena_model_t *bblgr2;

// modena_model_t *rheologymodel = modena_model_new("viscosity_Arrhenius");

modena_model_t *kinetics;

// modena_model_t *density_reaction_mixturemodel = modena_model_new("density_reaction_mixtureSM");
// modena_model_t *rheologymodel = modena_model_new("rheology");
// modena_model_t *kinetics = modena_model_new("simpleKinetics");
/* allocate memory and fetch arg positions:
    - bubbleGrowth1,
    - bubbleGrowth2,
    - density_reaction_mixture,
    - rheology.
*/
modena_inputs_t *inputs_bblgr1;
modena_outputs_t *outputs_bblgr1;

modena_inputs_t *inputs_bblgr2;
modena_outputs_t *outputs_bblgr2;

// modena_inputs_t *inputs_rheo       = modena_inputs_new (rheologymodel);
// modena_outputs_t *outputs_rheo     = modena_outputs_new (rheologymodel);

modena_inputs_t *inputs_kinetics;
modena_outputs_t *outputs_kinetics;

// modena_inputs_t *inputs_den        = modena_inputs_new (density_reaction_mixturemodel);
// modena_outputs_t *outputs_den      = modena_outputs_new (density_reaction_mixturemodel);

// modena_inputs_t *inputs_rheo       = modena_inputs_new (rheologymodel);
// modena_outputs_t *outputs_rheo     = modena_outputs_new (rheologymodel);

// modena_inputs_t *inputs_kinetics   = modena_inputs_new (kinetics);
// modena_outputs_t *outputs_kinetics = modena_outputs_new (kinetics);
#endif
