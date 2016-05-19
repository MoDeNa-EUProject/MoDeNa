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

size_t Catalyst_1_Pos;
size_t CE_A0_Pos;
size_t CE_A1_Pos;
size_t CE_B_Pos;
size_t CE_B2_Pos;
size_t CE_I0_Pos;
size_t CE_I1_Pos;
size_t CE_I2_Pos;
size_t CE_PBA_Pos;
size_t CE_Breac_Pos;
size_t CE_Areac0_Pos;
size_t CE_Areac1_Pos;
size_t CE_Ireac0_Pos;
size_t CE_Ireac1_Pos;
size_t CE_Ireac2_Pos;
size_t Bulk_Pos;
size_t R_1_Pos;
size_t R_1_mass_Pos;
size_t R_1_temp_Pos;
size_t R_1_vol_Pos;

// First 20 outputs (source terms)
size_t source_Catalyst_1_Pos;
size_t source_CE_A0_Pos;
size_t source_CE_A1_Pos;
size_t source_CE_B_Pos;
size_t source_CE_B2_Pos;
size_t source_CE_I0_Pos;
size_t source_CE_I1_Pos;
size_t source_CE_I2_Pos;
size_t source_CE_PBA_Pos;
size_t source_CE_Breac_Pos;
size_t source_CE_Areac0_Pos;
size_t source_CE_Areac1_Pos;
size_t source_CE_Ireac0_Pos;
size_t source_CE_Ireac1_Pos;
size_t source_CE_Ireac2_Pos;
size_t source_Bulk_Pos;
size_t source_R_1_Pos;
size_t source_R_1_mass_Pos;
size_t source_R_1_temp_Pos;
size_t source_R_1_vol_Pos;


#endif
