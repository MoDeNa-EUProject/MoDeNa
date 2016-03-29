/** @file modenaCalls.h
    @brief instantiates the surrogate models
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
if (modena_error_occurred())
{
    return modena_error();
}

bblgr2 = modena_model_new("bubbleGrowth2");
if (modena_error_occurred())
{
    return modena_error();
}

kinetics = modena_model_new("RF-1-public");
if (modena_error_occurred())
{
    return modena_error();
}

density_reaction_mixturemodel = modena_model_new("density_reaction_mixture");
if (modena_error_occurred())
{
    return modena_error();
}

inputs_bblgr1 = modena_inputs_new(bblgr1);
outputs_bblgr1 = modena_outputs_new(bblgr1);

inputs_bblgr2 = modena_inputs_new(bblgr2);
outputs_bblgr2 = modena_outputs_new(bblgr2);

inputs_kinetics = modena_inputs_new(kinetics);
outputs_kinetics = modena_outputs_new(kinetics);

inputs_den = modena_inputs_new (density_reaction_mixturemodel);
outputs_den = modena_outputs_new (density_reaction_mixturemodel);

#endif