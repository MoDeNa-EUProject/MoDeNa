/** @file readParameters.h
	@brief reads the inputs from the input files.
*/
void readParams();
void readParams() {
    std::ifstream inputs;
    inputs.open("../inputsQmom.in", std::ios::in);
    inputs >> Pr;   // initial/final pressure of the mixture, Pa
	inputs >> Temp0;		// initial temperature, K
/* Inputs for gelling reaction, XOH */
	inputs >> A_OH;		// m3/mol s
	inputs >> E_OH; 	// J/mol
	inputs >> OH_0;		// Initial concentration of polyol OH groups in the mixutre, mol/m3
	inputs >> NCO_0;		// Initial concentration of isocianate NCO groups in the mixutre, mol/m3
	inputs >> W_0;		// Initial concentration of water in the mixture, mol/m3
/* Inputs for blowing reaction, XW */
	inputs >> A_W;	// 1/s
	inputs >> E_W;	// J/mol
/* Constants */
	inputs >> RR;	// J/mol K
	inputs >> rhoPoly;	// kg/m3 Density of the liquid polymer
	inputs >> rhoBL;	// kg/m3 Density of the blowing agent
/* Inputs for enthalpy */
	inputs >> DH_OH;	// Reaction heat for the gelling reaction, J/mol
	inputs >> DH_W;	// Reaction heat for the blowing reaction, J/mol
	inputs >> C_Poly;	// Polyurethane specific heat, J/kg K
	inputs >> C_CO2;		// CO2 specific heat, J/kg K
	inputs >> C_BG;		// Physical blowing agent in gas phase specific heat, J/kg K
	inputs >> C_BL;		// Physical blowing agent in liquid phase specific heat, J/kg K
	inputs >> lambda;		// Latent heat of blowing agent, J/kg
// physical blowing agent (used for solubility model)
	inputs >>    phBL;			// 1=pentane, 2=R-11
// density model used
	inputs >>    denMod;			// 1=modena, 2=rhoPoly
// kinetics model used
	inputs >>    kinMod;			// 1=Baser, 2=Baser with R(x)
	inputs >>   dilution;		// use dilution effect
/* Inputs for weight fraction of gaseous CO2 in the mixture */
	inputs >> M_CO2;		// Molecular mass of carbon dioxide, kg/kmol
	inputs >> M_B;		// Molecular mass of blowing agent, kg/kmol
	inputs >> M_NCO;		// Molecular weight of NCO, kg/kmol
	inputs >> M_air;		// Molecular weight of air, kg/kmol
	inputs >> CO2_D;	// Weight fraction of dissolved CO2 in the mixture, -
	inputs >> L0;	// Initial weight fraction of blowing agent in the liquid, -
	inputs >> CO2_0;		// Initial weight fraction of CO2 in the liquid, -
// Other physical properties
	inputs >> surfaceTension; // required for partial pressure
// initial bubble size distribution
	inputs >> sig;		// correlated to variance of initial distribution
	inputs >> init_size;	// initial mean bubble diameter, m
	inputs >> NN;	// correlated to number of initial bubbles in m^3
// integration parameters
	inputs >> abs_err;// absolute error
	inputs >> rel_err;// relative error
	inputs >> dt;		// time step, s
	inputs >> tend;		// end time, s
    inputs.close();
    if (W_0<1e-8) {
        W_0=-1.0;
    }
}
