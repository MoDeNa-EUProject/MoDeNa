/**
* @ingroup mod_3Dcfd
* @file modenaData.H
 @brief allocate memory for the surrogate models
*/

#ifndef MODENADATA_H
#define MODENADATA_H
/// pointer to first surrogate model for the bubble growth
modena_model_t *bblgr1;
/// pointer to second surrogate model for the bubble growth
modena_model_t *bblgr2;
/// pointer to the surrogate model for kinetics of reactions
modena_model_t *kinetics;
/// pointer to the surrogate model for density of reacting mixture
modena_model_t *density_reaction_mixturemodel;
/// pointer to the surrogate model for rheology model
modena_model_t *rheologymodel;
/// pointer to the surrogate model for strut content
modena_model_t *strutContentmodel;
/// pointer to the surrogate model for thermal conductivity
modena_model_t *thermalConductivitymodel;
/// pointer for the inputs of first bubble growth surrogate model
modena_inputs_t *inputs_bblgr1;
/// pointer for the outputs of first bubble growth surrogate model
modena_outputs_t *outputs_bblgr1;
/// pointer for the inputs of second bubble growth surrogate model
modena_inputs_t *inputs_bblgr2;
/// pointer for the outputs of second bubble growth surrogate model
modena_outputs_t *outputs_bblgr2;
/// pointer for the inputs of kinetics surrogate model
modena_inputs_t *inputs_kinetics;
/// pointer for the outputs of kinetics surrogate model
modena_outputs_t *outputs_kinetics;
/// pointer for the inputs of density of reacting mixture surrogate model
modena_inputs_t *inputs_den;
/// pointer for the outputs of density of reacting mixture surrogate model
modena_outputs_t *outputs_den;
/// pointer for the inputs of rheology surrogate model
modena_inputs_t *inputs_rheo;
/// pointer for the outputs of rheology surrogate model
modena_outputs_t *outputs_rheo;
/// pointer for the inputs of strut content surrogate model
modena_inputs_t *inputs_strutContent;
/// pointer for the outputs of strut content surrogate model
modena_outputs_t *outputs_strutContent;
/// pointer for the inputs of thermal conductivity surrogate model
modena_inputs_t *inputs_thermalConductivity;
/// pointer for the outputs of thermal conductivity surrogate model
modena_outputs_t *outputs_thermalConductivity;


/// @param[in] size_t kineticTime_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t kineticTime_Pos;
/// @param[in] size_t Catalyst_1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t Catalyst_1_Pos;
/// @param[in] size_t CE_A0_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_A0_Pos;
/// @param[in] size_t CE_A1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_A1_Pos;
/// @param[in] size_t CE_B_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_B_Pos;
/// @param[in] size_t CE_B2_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_B2_Pos;
/// @param[in] size_t CE_I0_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_I0_Pos;
/// @param[in] size_t CE_I1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_I1_Pos;
/// @param[in] size_t CE_I2_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_I2_Pos;
/// @param[in] size_t CE_PBA_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_PBA_Pos;
/// @param[in] size_t CE_Breac_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_Breac_Pos;
/// @param[in] size_t CE_Areac0_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_Areac0_Pos;
/// @param[in] size_t CE_Areac1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_Areac1_Pos;
/// @param[in] size_t CE_Ireac0_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_Ireac0_Pos;
/// @param[in] size_t CE_Ireac1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_Ireac1_Pos;
/// @param[in] size_t CE_Ireac2_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t CE_Ireac2_Pos;
/// @param[in] size_t Bulk_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t Bulk_Pos;
/// @param[in] size_t R_1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t R_1_Pos;
/// @param[in] size_t R_1_temp_RF1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t R_1_temp_RF1_Pos;
/// @param[in] size_t R_1_mass_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t R_1_mass_Pos;
/// @param[in] size_t R_1_vol_RF1_Pos
/// @brief memory allocation for inputs of kinetics surrogate model
size_t R_1_vol_RF1_Pos;
/// @param[out] size_t source_Catalyst_1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_Catalyst_1_Pos;
/// @param[out] size_t source_CE_A0_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_A0_Pos;
/// @param[out] size_t source_CE_A1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_A1_Pos;
/// @param[out] size_t source_CE_B_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_B_Pos;
/// @param[out] size_t source_CE_B2_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_B2_Pos;
/// @param[out] size_t source_CE_I0_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_I0_Pos;
/// @param[out] size_t source_CE_I1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_I1_Pos;
/// @param[out] size_t source_CE_I2_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_I2_Pos;
/// @param[out] size_t source_CE_PBA_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_PBA_Pos;
/// @param[out] size_t source_CE_Breac_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_Breac_Pos;
/// @param[out] size_t source_CE_Areac0_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_Areac0_Pos;
/// @param[out] size_t source_CE_Areac1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_Areac1_Pos;
/// @param[out] size_t source_CE_Ireac0_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_Ireac0_Pos;
/// @param[out] size_t source_CE_Ireac1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_Ireac1_Pos;
/// @param[out] size_t source_CE_Ireac2_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_CE_Ireac2_Pos;
/// @param[out] size_t source_Bulk_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_Bulk_Pos;
/// @param[out] size_t source_R_1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_R_1_Pos;
/// @param[out] size_t source_R_1_mass_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_R_1_mass_Pos;
/// @param[out] size_t source_R_1_temp_RF1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_R_1_temp_RF1_Pos;
/// @param[out] size_t source_R_1_vol_RF1_Pos
/// @brief memory allocation for outputs of kinetics surrogate model
size_t source_R_1_vol_RF1_Pos;

/// @param[in] size_t rho_foam_Pos
/// @brief memory allocation for the foam density as the input of thermal conductivity surrogate model
size_t rho_foam_Pos;
/// @param[in] size_t porosity_Pos
/// @brief memory allocation for the gas volume fraction as the input of thermal conductivity surrogate model
size_t porosity_Pos;
/// @param[in] size_t cell_size_Pos
/// @brief memory allocation for the bubble size as the input of thermal conductivity surrogate model
size_t cell_size_Pos;
/// @param[in] size_t strut_c_Pos
/// @brief memory allocation for the strut content as the input of thermal conductivity surrogate model
size_t strut_c_Pos;
/// @param[in] size_t temp_Pos
/// @brief memory allocation for the temperature as the input of thermal conductivity surrogate model
size_t temp_Pos;
/// @param[in] size_t X_CO2_Pos
/// @brief memory allocation for the partial pressure of CO2 as the input of thermal conductivity surrogate model
size_t X_CO2_Pos;
/// @param[in] size_t X_Cyp_Pos
/// @brief memory allocation for the partial pressure of cyclo-pentane as the input of thermal conductivity surrogate model
size_t X_Cyp_Pos;
/// @param[in] size_t X_O2_Pos
/// @brief memory allocation for the partial pressure of O2 as the input of thermal conductivity surrogate model
size_t X_O2_Pos;
/// @param[in] size_t X_N2_Pos
/// @brief memory allocation for the partial pressure of N2 as the input of thermal conductivity surrogate model
size_t X_N2_Pos;
/// @param[in] size_t strut_content_Pos
/// @brief memory allocation for the strut content as the input of thermal conductivity surrogate model
size_t strut_content_Pos;

// rheology MoDeNa
/// @param[in] size_t temp_rheopos
/// @brief memory allocation for the temperature as the input of rheology surrogate model
size_t temp_rheopos;
/// @param[in] size_t shear_rheopos
/// @brief memory allocation for the shear rate as the input of rheology surrogate model
size_t shear_rheopos;
/// @param[in] size_t conv_rheopos
/// @brief memory allocation for the conversion as the input of rheology surrogate model
size_t conv_rheopos;
/// @param[in] size_t m0_rheopos
/// @brief memory allocation for the moment of order zero as the input of rheology surrogate model
size_t m0_rheopos;
/// @param[in] size_t m1_rheopos
/// @brief memory allocation for the moment of order one as the input of rheology surrogate model
size_t m1_rheopos;
#endif
