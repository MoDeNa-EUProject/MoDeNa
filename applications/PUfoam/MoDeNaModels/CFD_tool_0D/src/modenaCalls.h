/** @file modenaCalls.h
	@brief instantiates the surrogate models, allocates memory and fetches arg positions
*/
#ifndef MODENACALLS_H
#define MODENACALLS_H
/*
instantiate the surrogate models:
    - bubbleGrowth1,
    - bubbleGrowth2,
    - density_reaction_mixture,
    - rheology,
    - simpleKinetics
*/

bblgr1 = modena_model_new("bubbleGrowth1");
if(modena_error_occurred())
{
    return modena_error();
}

bblgr2 = modena_model_new("bubbleGrowth2");
if(modena_error_occurred())
{
    return modena_error();
}

// modena_model_t *rheologymodel = modena_model_new("viscosity_Arrhenius");

kinetics = modena_model_new("RF-1-public");
if(modena_error_occurred())
{
    return modena_error();
}

// modena_model_t *density_reaction_mixturemodel = modena_model_new("density_reaction_mixtureSM");
// modena_model_t *rheologymodel = modena_model_new("rheology");
// modena_model_t *kinetics = modena_model_new("simpleKinetics");
/* allocate memory and fetch arg positions:
    - bubbleGrowth1,
    - bubbleGrowth2,
    - density_reaction_mixture,
    - rheology.
*/
inputs_bblgr1     = modena_inputs_new(bblgr1);
outputs_bblgr1   = modena_outputs_new(bblgr1);

inputs_bblgr2     = modena_inputs_new(bblgr2);
outputs_bblgr2   = modena_outputs_new(bblgr2);

// modena_inputs_t *inputs_rheo       = modena_inputs_new (rheologymodel);
// modena_outputs_t *outputs_rheo     = modena_outputs_new (rheologymodel);

inputs_kinetics   = modena_inputs_new(kinetics);
outputs_kinetics = modena_outputs_new(kinetics);

// modena_inputs_t *inputs_den        = modena_inputs_new (density_reaction_mixturemodel);
// modena_outputs_t *outputs_den      = modena_outputs_new (density_reaction_mixturemodel);

// modena_inputs_t *inputs_rheo       = modena_inputs_new (rheologymodel);
// modena_outputs_t *outputs_rheo     = modena_outputs_new (rheologymodel);

// modena_inputs_t *inputs_kinetics   = modena_inputs_new (kinetics);
// modena_outputs_t *outputs_kinetics = modena_outputs_new (kinetics);
#endif
