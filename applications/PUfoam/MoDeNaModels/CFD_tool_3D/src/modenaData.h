/** @file modenaCalls.h
	@brief Allocate memory for the surrogate models
*/
#ifndef MODENADATA_H
#define MODENADATA_H

modena_model_t *bblgr1;
modena_model_t *bblgr2;

modena_model_t *kinetics;

modena_model_t *density_reaction_mixturemodel;

modena_model_t *rheologymodel;

modena_model_t *strutContentmodel;

modena_model_t *thermalConductivitymodel;

modena_inputs_t *inputs_bblgr1;
modena_outputs_t *outputs_bblgr1;

modena_inputs_t *inputs_bblgr2;
modena_outputs_t *outputs_bblgr2;

modena_inputs_t *inputs_kinetics;
modena_outputs_t *outputs_kinetics;

modena_inputs_t *inputs_den;
modena_outputs_t *outputs_den;

modena_inputs_t *inputs_rheo;
modena_outputs_t *outputs_rheo;

modena_inputs_t *inputs_strutContent;
modena_outputs_t *outputs_strutContent;

modena_inputs_t *inputs_thermalConductivity;
modena_outputs_t *outputs_thermalConductivity;

size_t kineticTime_Pos;
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
size_t R_1_temp_RF1_Pos;
size_t R_1_mass_Pos;
size_t R_1_vol_RF1_Pos;
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
size_t source_R_1_temp_RF1_Pos;
size_t source_R_1_vol_RF1_Pos;

// strut content
size_t rho_foam_Pos;
size_t strut_content_Pos;

// thermal conductivity
size_t porosity_Pos;
size_t cell_size_Pos;
size_t strut_c_Pos;
size_t temp_Pos;
size_t X_CO2_Pos;
size_t X_Cyp_Pos;
size_t X_O2_Pos;
size_t X_N2_Pos;

// rheology MoDeNa
size_t temp_rheopos;
size_t shear_rheopos;
size_t conv_rheopos;
size_t m0_rheopos;
size_t m1_rheopos;
#endif
