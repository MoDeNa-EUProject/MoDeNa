/* Description:
	All the experimental data required for the source terms calculations
*/
	double Pr		= 101325.0; // initial/final pressure of the mixture, Pa
	double Temp0	= 297;		// initial temperature, K
/* Inputs for gelling reaction, XOH */
	double A_OH     = 1965.0;		// m3/mol s
	double E_OH     = 5.487e4; 	// J/mol
	double OH_0   = 3765.0;		// Initial concentration of polyol OH groups in the mixutre, mol/m3
	double NCO_0  = 3765.0;		// Initial concentration of isocianate NCO groups in the mixutre, mol/m3
	double W_0    = -1.0;		// Initial concentration of water in the mixture, mol/m3
/* Inputs for blowing reaction, XW */
	double A_W      = 1.385e3;	// 1/s
	double E_W      = 3.266e4;	// J/mol
/* Constants */
	double RR      = 8.3145;	// J/mol K
	double rhoPoly = 1100.0;	// kg/m3 Density of the liquid polymer
	double rhoBL   = 1467.0;	// kg/m3 Density of the blowing agent
/* Inputs for enthalpy */
	double DH_OH   = -7.49e4;	// Reaction heat for the gelling reaction, J/mol
	double DH_W    = -8.6e4;	// Reaction heat for the blowing reaction, J/mol
	double C_Poly  = 1800.0;	// Polyurethane specific heat, J/kg K
	double C_CO2   = 836.6;		// CO2 specific heat, J/kg K
	double C_BG    = 593.0;		// Physical blowing agent in gas phase specific heat, J/kg K
	double C_BL    = 870.0;		// Physical blowing agent in liquid phase specific heat, J/kg K
	double lambda  = 2.0e5;		// Latent heat of blowing agent, J/kg
	double C_TOT;				// Total specifc heat of the mixture
// physical blowing agent (used for solubility model)
	int    phBL    = 2;			// 1=pentane, 2=R-11
// density model used
	int    denMod  = 2;			// 1=modena, 2=rhoPoly
// kinetics model used
	int    kinMod  = 2;			// 1=Baser, 2=Baser with R(x)
	bool   dilution= true;		// use dilution effect
/* Inputs for weight fraction of gaseous CO2 in the mixture */
	double M_CO2   = 44.0;		// Molecular mass of carbon dioxide, kg/kmol
	double M_B     = 137.37;		// Molecular mass of blowing agent, kg/kmol
	double M_NCO   = 615.0;		// Molecular weight of NCO, kg/kmol
	double M_air   = 29.0;		// Molecular weight of air, kg/kmol
	double CO2_D   = 4.4e-4;	// Weight fraction of dissolved CO2 in the mixture, -
	double L0      = 0.155;	// Initial weight fraction of blowing agent in the liquid, -
	double CO2_0   = 0.0;		// Initial weight fraction of CO2 in the liquid, -
// Other physical properties
	double surfaceTension = 25e-3; // required for partial pressure
// initial bubble size distribution
	double sig 		= 1e-1;		// correlated to variance of initial distribution
	double init_size= 2.0e-6;	// initial mean bubble diameter, m
	double NN 		= 0.138e12;	// correlated to number of initial bubbles in m^3
	double air_g;				// air weight fraction
// integration parameters
	double abs_err = 1e-6;// absolute error
	double rel_err = 1e-6;// relative error
	double dt = 1e0;		// time step, s
	double tend = 200.0;		// end time, s
