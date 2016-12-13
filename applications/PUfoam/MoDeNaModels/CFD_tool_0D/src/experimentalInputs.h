/** @file experimentalInputs.h
	@brief All the experimental data required for the source terms calculations
	@var double Pr
	@brief initial/final pressure of the mixture, Pa
	@var double Temp0
	@brief initial temperature, K
	@var double A_OH
	@brief pre-exponential factor for the gelling reaction, 1/s
	@var double E_OH
	@brief activation energy for the gelling reaction, J/mol
	@var double OH_0
	@brief Initial concentration of polyol OH groups in the mixutre, mol/m3
	@var double NCO_0
	@brief Initial concentration of isocianate NCO groups in the mixutre, mol/m3
	@var double W_0
	@brief Initial concentration of water in the mixture, mol/m3
	@var double A_W
	@brief pre-exponential factor for the blowing reaction, 1/s
	@var double E_W
	@brief activation energy for the blowing reaction, J/mol
	@var double RR
	@brief ideal gas constant, J/mol K
	@var double rhoPoly
	@brief density of the liquid polymer, kg/m3
	@var double rhoBL
	@brief density of the blowing agent, kg/m3
	@var double DH_OH
	@brief Reaction heat for the gelling reaction, J/mol
	@var double DH_W
	@brief Reaction heat for the blowing reaction, J/mol
	@var double C_Poy
	@brief Polyurethane specific heat, J/kg K
	@var double C_CO2
	@brief CO2 specific heat, J/kg K
	@var double C_BG
	@brief Specific heat of the physical blowing agent in gas phase, J/kg K
	@var double C_BL
	@brief Specific heat of the physical blowing agent in liquid phase, J/kg K
	@var double lambda
	@brief Latent heat of blowing agent, J/kg
	@var double C_TOT
	@brief Total specifc heat of the mixture
	@var int phBL
	@brief type of physical blowing agent 1 = pentane, 2 = R-11
	@var int denMod
	@brief density mode, 1 = modena, 2 = constant
	@var int kinMod
	@brief kinetics model, 1 = Baser, 2 = Baser with R(x)
	@sa http://onlinelibrary.wiley.com/doi/10.1002/pen.760340804/abstract
	@var bool dilution
	@brief use dilution effect
	@var double M_CO2
	@brief Molecular mass of carbon dioxide, kg/kmol
	@var double M_B
	@brief Molecular mass of blowing agent, kg/kmol
	@var double M_NCO
	@brief Molecular weight of NCO, kg/kmol
	@var double M_air
	@brief Molecular weight of air, kg/kmol
	@var double CO2_D
	@brief Weight fraction of dissolved CO2 in the mixture, -
	@var double L0
	@brief Initial weight fraction of blowing agent in the liquid, -
	@var double CO2_0
	@brief Initial weight fraction of CO2 in the liquid, -
	@var double surfaceTension
	@brief required for the computation of partial pressure
	@var double sig
	@brief correlated to variance of initial distribution
	@var double init_size
	@brief initial mean bubble diameter, m
	@var double NN
	@brief correlated to number of initial bubbles in m^3
	@var double air_g
	@brief air weight fraction
	@var double abs_err
	@brief absolute error
	@var double rel_err
	@brief relative error
	@var double dt
	@brief time step, s
	@var double tend
	@brief  end time, s
*/
	double Pr		= 101325.0;
	double Temp0	= 297;
/* Inputs for gelling reaction, XOH */
	double A_OH     = 1965.0;
	double E_OH     = 5.487e4;
	double OH_0   = 3765.0;
	double NCO_0  = 3765.0;
	double W_0    = -1.0;

	double catalyst;
	double polyol1_ini;
	double polyol2_ini;
	double amine_ini;
	double isocyanate1_ini;
	double isocyanate2_ini;
	double isocyanate3_ini;

/* Inputs for blowing reaction, XW */
	double A_W      = 1.385e3;
	double E_W      = 3.266e4;
/* Constants */
	double RR      = 8.3145;
	double rhoPoly = 1100.0;
	double rhoBL   = 1467.0;
/* Inputs for enthalpy */
	double DH_OH   = -7.49e4;
	double DH_W    = -8.6e4;
	double C_Poly  = 1800.0;
	double C_CO2   = 836.6;
	double C_BG    = 593.0;
	double C_BL    = 870.0;
	double lambda  = 2.0e5;
	double C_TOT;
// physical blowing agent (used for solubility model)
	int    phBL    = 1;		// 1=pentane, 2=R-11
// density model used
	int    denMod;
// kinetics model used
	int    kinMod  = 2;
	bool   dilution= true;
	double X_gel;
/* Inputs for weight fraction of gaseous CO2 in the mixture */
	double M_CO2   = 44.0;
	double M_B     = 137.37;
	double M_NCO   = 615.0;
	double M_air   = 29.0;
	double L0      = 0.155;
	double CO2_0   = 0.0;
// Other physical properties
	double surfaceTension = 25e-3;
// initial bubble size distribution
	double sig 		= 1e-1;
	double init_size= 2.0e-6;
	double NN 		= 0.138e12;
	double air_g;
// integration parameters
	double abs_err = 1e-6;
	double rel_err = 1e-6;
	double dt = 1e0;
	double tend = 200;
// Realizability
	bool realizabilityCheck = false; 	// switch for realizability
// 2nodes vs meanDiameter
	std::string bubbleMode; 	// mean radius, two nodes
	bool apparentViscosity;
	// bool kinetics_basf = false;
