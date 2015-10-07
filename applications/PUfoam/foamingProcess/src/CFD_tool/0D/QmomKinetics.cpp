/**
@cond
   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2015 MoDeNa Consortium, All rights reserved.

License
    This file is part of Modena.

    Modena is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.
@endcond
@file
    This is a macro-scale modeling tool for the foaming process. The code utilizes 
    the MoDeNa interface library to connect different models including nano, and 
    meso scale models. The code returns the evolution of foam properties such as
    density, temperature and bubble/cell size distribution.
@brief macor-scale tool for the foaming process.   
@authors    Mohsen Karimi, Daniele Marchisio, Pavel Ferkl
@copyright  2014-2015, MoDeNa Project. GNU Public License.
@ingroup    app_foaming
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "modena.h"

extern "C"{void dsteqr_(char &, int *, double *, double *, double *, int *, double *, int *); }

#include "experimentalInputs.h"
#include "readParameters.h"
#include "initializeMoments.h"
#include "pda.h"
#include "growth.h"
#include "coalescence.h"
#include "liquidBA.h"

using namespace std;
using namespace boost::numeric::odeint;
/**
@typedef 
typedef vector<double> to state_type 
typedef runge_kutta_cash_karp54< state_type > error_stepper_type
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type
*/
typedef std::vector< double > state_type;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
/**
@var controlled_stepper
@sa http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/concepts/controlled_stepper.html 
*/
controlled_stepper_type controlled_stepper;

#include "partialPressure.h"
/**
@var dpdt[2]: global double array 
@brief This is used to compute the partial pressures.

@var pOld[2]: global double array variable 
@brief This is to hold the old pressure values during the partial pressure calculations.
*/
double dpdt[2] = {};
double pOld[2] = {};

#include "modenaCalls.h"
#include "momentsConverter.h"
#include "write_kinetics.h"
/**
@fn QmomKinetics(const state_type &y , state_type &dydt , double t)
@brief This is to calculate the RHD of all the ODEs.
@param [in] const state_type &y - vector<double> to hold the results.
@param [in] state_type &dydt -  vector<double> to hold the RHDs of ODEs.
@param [in] double t - time
@return void.
*/
void QmomKinetics( const state_type &y , state_type &dydt , double t )
{
	// dydt[0] : XW
	// dydt[1] : XOH
	// dydt[2] : T
	// dydt[3] : L_l
	// dydt[4] : L_g
	// dydt[5] : CO2_l
	// dydt[6] : CO2_g
	// dydt[7] : m0
	// dydt[8] : m1
	// dydt[9] : m2
	// dydt[10]: m3
	// ---simpleKinetics---
    // dydt[11]: EG_NCO
    // dydt[12]: EG_OH
    // dydt[13]: H2O
    // dydt[14]: CO2
    // dydt[15]: PENTANE
    // dydt[16]: POLYMER
    // dydt[17]: POLYMERBLOW
    // dydt[18]: UREA
    // dydt[19]: R_1_temp

	int nNodes = 2;
	int mOrder[6] = {0, 1, 2, 3, 4, 5};
	double mom[2*nNodes], we[nNodes], vi[nNodes], sgBA[2*nNodes], sgCO2[2*nNodes], sc[2*nNodes];
	double L_l, L_g, CO2_l, CO2_g, T, Lm, c1;
	double rhoPolySurrgate;
	double beta0 	= 0.0;
    double XW, XOH;
	// simpleKinetics variables
    double EG_XNCO, EG_NCO, EG_OH, EG_XOH, H2O, XH2O, CO2, PENTANE, POLYMER, POLYMERBLOW, UREA, R_1_temp;
    double init_EG_NCO 	= 5.0;
    double init_EG_OH 	= 5.0;
    double init_H2O 	= 0.2;

    XW      	= y[0];
    XOH     	= y[1];
    T       	= y[2];
    L_l     	= y[3];
    L_g     	= y[4];
    CO2_l   	= y[5];
    CO2_g   	= y[6];
    mom[0]  	= y[7];
    mom[1]  	= y[8];
    mom[2]  	= y[9];
    mom[3]  	= y[10];
    // EG_NCO 		= y[11];
    // EG_OH 		= y[12];
    // H2O 		= y[13];
    // CO2         = y[14];
    // PENTANE     = y[15];
    // POLYMER     = y[16];
    // POLYMERBLOW = y[17];
    // UREA        = y[18];
    // R_1_temp    = y[19];

    // EG_XNCO     = 1.0 - (y[11]/init_EG_NCO);
    // EG_XOH      = 1.0 - (y[12]/init_EG_OH);
    // XH2O        = 1.0 - (y[13]/init_H2O);

    // Calling the simpleKinetics model
	// if (kinMod == 3)
 //    {
	//     // inputs argPos
	// 	modena_model_t *kinetics = modena_model_new("simpleKinetics");
	// 	modena_inputs_t *inputs_kinetics   = modena_inputs_new (kinetics);
	// 	modena_outputs_t *outputs_kinetics = modena_outputs_new (kinetics);
	//     size_t EG_NCO_Pos       = modena_model_inputs_argPos(kinetics, "'EG_NCO'");
	//     size_t EG_OH_Pos        = modena_model_inputs_argPos(kinetics, "'EG_OH'");
	//     size_t H2O_Pos          = modena_model_inputs_argPos(kinetics, "'H2O'");
	//     size_t CO2_Pos          = modena_model_inputs_argPos(kinetics, "'CO2'");
	//     size_t PENTANE_Pos      = modena_model_inputs_argPos(kinetics, "'PENTANE'");
	//     size_t POLYMER_Pos      = modena_model_inputs_argPos(kinetics, "'POLYMER'");
	//     size_t POLYMERBLOW_Pos  = modena_model_inputs_argPos(kinetics, "'POLMERBLOW'");
	//     size_t UREA_Pos         = modena_model_inputs_argPos(kinetics, "'UREA'");
	//     size_t R_1_temp_Pos     = modena_model_inputs_argPos(kinetics, "'R_1_temp'");

	//     // outputs argPos
	//     size_t source_EG_NCO_Pos        = modena_model_outputs_argPos(kinetics, "source_EG_NCO");
	//     size_t source_EG_OH_Pos         = modena_model_outputs_argPos(kinetics, "source_EG_OH");
	//     size_t source_H2O_Pos           = modena_model_outputs_argPos(kinetics, "source_H2O");
	//     size_t source_CO2_Pos           = modena_model_outputs_argPos(kinetics, "source_CO2");
	//     size_t source_PENTANE_Pos       = modena_model_outputs_argPos(kinetics, "source_PENTANE");
	//     size_t source_POLYMER_Pos       = modena_model_outputs_argPos(kinetics, "source_POLYMER");
	//     size_t source_POLYMERBLOW_Pos   = modena_model_outputs_argPos(kinetics, "source_POLMERBLOW");
	//     size_t source_UREA_Pos          = modena_model_outputs_argPos(kinetics, "source_UREA");
	//     size_t source_R_1_temp_Pos      = modena_model_outputs_argPos(kinetics, "source_R_1_temp");

	//     modena_model_argPos_check(kinetics);

	//     // set input vector
	//     modena_inputs_set(inputs_kinetics, EG_NCO_Pos, EG_NCO);
	//     modena_inputs_set(inputs_kinetics, EG_OH_Pos, EG_OH);
	//     modena_inputs_set(inputs_kinetics, H2O_Pos, H2O);
	//     modena_inputs_set(inputs_kinetics, CO2_Pos, CO2);
	//     modena_inputs_set(inputs_kinetics, PENTANE_Pos, PENTANE);
	//     modena_inputs_set(inputs_kinetics, POLYMER_Pos, POLYMER);
	//     modena_inputs_set(inputs_kinetics, POLYMERBLOW_Pos, POLYMERBLOW);
	//     modena_inputs_set(inputs_kinetics, UREA_Pos, UREA);
	//     modena_inputs_set(inputs_kinetics, R_1_temp_Pos, R_1_temp);


	//     // call the model
	//     int ret_kinetics = modena_model_call(kinetics, inputs_kinetics, outputs_kinetics);

	//     // terminate, if requested
	//     if(ret_kinetics != 0)
	//     {
	//         modena_inputs_destroy (inputs_kinetics);
	//         modena_outputs_destroy (outputs_kinetics);
	//         modena_model_destroy (kinetics);
	//         //return ret_kinetics;
	//     }

	// 	// get the source terms for simpleKinetics
	//     dydt[11] = modena_outputs_get(outputs_kinetics, source_EG_NCO_Pos);
	//     dydt[12] = modena_outputs_get(outputs_kinetics, source_EG_OH_Pos);
	//     dydt[13] = modena_outputs_get(outputs_kinetics, source_H2O_Pos);
	//     dydt[14] = modena_outputs_get(outputs_kinetics, source_CO2_Pos);
	//     dydt[15] = modena_outputs_get(outputs_kinetics, source_PENTANE_Pos);
	//     dydt[16] = modena_outputs_get(outputs_kinetics, source_POLYMER_Pos);
	//     dydt[17] = modena_outputs_get(outputs_kinetics, source_POLYMERBLOW_Pos);
	//     dydt[18] = modena_outputs_get(outputs_kinetics, source_UREA_Pos);
	//     dydt[19] = modena_outputs_get(outputs_kinetics, source_R_1_temp_Pos);
	// }

	// Check for negative sources
	// for (int i = 11; i < 20; i++)
	// {
	// 	if(dydt[i] < 0.0)
	// 	{
	// 		dydt[i] = 0.0;
	// 	}
	// }

	switch (denMod)
    {
		case 1:
        {
            // Calling the model for density reaction mixture
            // size_t T_denpos     = modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
            // size_t XOH_denpos   = modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
            // modena_model_argPos_check(density_reaction_mixturemodel);

            // // set input vector
            // modena_inputs_set(inputs_den, T_denpos, T);
            // modena_inputs_set(inputs_den, XOH_denpos, EG_XOH);

            // // call the model
            // int ret_den = modena_model_call (density_reaction_mixturemodel, inputs_den, outputs_den);
            
            // if(ret_den != 0)
            // {
            //     modena_inputs_destroy (inputs_den);
            //     modena_outputs_destroy (outputs_den);
            //     modena_model_destroy (density_reaction_mixturemodel);
            //     exit(ret_den);
            // }

            // rhoPolySurrgate = modena_outputs_get(outputs_den, 0);
            // break;
		}
		case 2:
			rhoPolySurrgate = rhoPoly;
            break;
	}

	// ODEs
    dydt[0] 	= A_W*exp(-E_W/(RR*y[2]))*(1-y[0]);
    if(dydt[0] < 0.0)
    {
    	dydt[0] = 0.0;
    }
    if(W_0 < 0.0)
    {
    	dydt[0] = 0.0;
    }

	double Rx;
	switch (kinMod) {
		case 1:
			Rx=1;
		case 2:
			if (y[1]<0.5) {
				Rx=1;
			} else if (y[1]<0.87) {
				Rx=-2.027*y[1]+2.013;
			} else {
				Rx=3.461*y[1]-2.761;
			}
	}

    dydt[1] 	= Rx*A_OH*exp(-E_OH/(RR*y[2]))*OH_0*(1-y[1])*(NCO_0/OH_0 - 2.0*y[0]*W_0/OH_0 - y[1]);
    if(dydt[1] < 0.0)
    {
    	dydt[1] = 0.0;
    }
	if (dilution) {
		dydt[0]=dydt[0]/(1+y[3]*rhoPolySurrgate/rhoBL);
		dydt[1]=dydt[1]/(1+y[3]*rhoPolySurrgate/rhoBL);
	}
	double dT=1e-4;
	double dLdT=(min(LMax(T),L0)-min(LMax(T+dT),L0))/dT;
    C_TOT 		= C_Poly + CO2_g*C_CO2 + L_g*C_BG + L_l*C_BL + dLdT*lambda;
		// this implementation of evaporation heat assumes that concentration of
		// physical blowing agent in liquid is always in equilibrium
    dydt[2] 	= (-DH_OH*OH_0)/(rhoPolySurrgate*C_TOT)*dydt[1]+(-DH_W*W_0)/(rhoPolySurrgate*C_TOT)*dydt[0];

	Lm 			= LMax(T);

	// bubble radius for bblgr1 and bblgr2 model
    double R 	= bubbleRadius(mom[0], mom[1]);

    // partial pressure within bubbles due to the evaporation of physical blowing agent
    double p_1  = partialPressureBA(y);

    // partial pressure within bubbles due to the generation of CO2
    double p_2  = partialPressureCO2(y);

    double c_1  = L_l*rhoPolySurrgate*1000.0/M_B;
    double c_2  = CO2_l*rhoPolySurrgate*1000.0/M_CO2;
    double KH1  = (rhoPolySurrgate*Lm)/((M_B/1000.0)*Pr);
    double KH2  = (rhoPolySurrgate*CO2_D)/((M_CO2/1000.0)*Pr);

    size_t Tbblgr1pos               = modena_model_inputs_argPos(bblgr1, "T");
    size_t Rbblgr1pos 	            = modena_model_inputs_argPos(bblgr1, "R");
    size_t KH1bblgr1pos 	        = modena_model_inputs_argPos(bblgr1, "kH");
    size_t c_1bblgr1pos             = modena_model_inputs_argPos(bblgr1, "c");
    size_t p_1bblgr1pos             = modena_model_inputs_argPos(bblgr1, "p");
    modena_model_argPos_check(bblgr1);
    size_t Tbblgr2pos               = modena_model_inputs_argPos(bblgr2, "T");
    size_t Rbblgr2pos               = modena_model_inputs_argPos(bblgr2, "R");
    size_t KH2bblgr2pos             = modena_model_inputs_argPos(bblgr2, "kH");
    size_t c_2bblgr2pos             = modena_model_inputs_argPos(bblgr2, "c");
    size_t p_2bblgr2pos             = modena_model_inputs_argPos(bblgr2, "p");
    modena_model_argPos_check(bblgr2);

    // set input vector
    modena_inputs_set(inputs_bblgr1, Tbblgr1pos, T);
    modena_inputs_set(inputs_bblgr1, Rbblgr1pos, R);
    modena_inputs_set(inputs_bblgr1, KH1bblgr1pos, KH1);
    modena_inputs_set(inputs_bblgr1, c_1bblgr1pos, c_1);
    modena_inputs_set(inputs_bblgr1, p_1bblgr1pos, p_1);
	// set input vector
    modena_inputs_set(inputs_bblgr2, Tbblgr2pos, T);
    modena_inputs_set(inputs_bblgr2, Rbblgr2pos, R);
    modena_inputs_set(inputs_bblgr2, KH2bblgr2pos, KH2);
    modena_inputs_set(inputs_bblgr2, c_2bblgr2pos, c_2);
    modena_inputs_set(inputs_bblgr2, p_2bblgr2pos, p_2);

    // call the bblgr1 model
    int ret_bblgr1 = modena_model_call (bblgr1, inputs_bblgr1, outputs_bblgr1);
    // terminate, if requested
    if(ret_bblgr1 != 0)
    {
        modena_inputs_destroy (inputs_bblgr1);
        modena_outputs_destroy (outputs_bblgr1);
        modena_model_destroy (bblgr1);
        // return ret_bblgr1;
    }
	// call the bblgr2 model
    int ret_bblgr2 = modena_model_call (bblgr2, inputs_bblgr2, outputs_bblgr2);
    // terminate, if requested
    if(ret_bblgr2 != 0)
    {
        modena_inputs_destroy (inputs_bblgr2);
        modena_outputs_destroy (outputs_bblgr2);
        modena_model_destroy (bblgr2);
        // return ret_bblgr2;
    }

    double G1, G2, dVdt_1, dVdt_2;
    G1 = modena_outputs_get(outputs_bblgr1, 0);
    G2 = modena_outputs_get(outputs_bblgr2, 0);
	// double mpar=0.0;
	// double mpar2=1.0;
	// G1=G1*pow(R,mpar)*mpar2; //for testing
	// G2=G2*pow(R,mpar)*mpar2;
	dVdt_1 = (G1*RR*T)/(p_1);
    if (dVdt_1 < 0.0 || G1 < 0.0 || L0<1e-8 || y[1]>0.5) //hardcoded gel point
    {
    	dVdt_1 = 0.0;
    }
	dVdt_2 = (G2*RR*T)/(p_2);
	if (dVdt_2 < 0.0 || G2 < 0.0 || W_0<1e-8 || y[1]>0.5) //hardcoded gel point
    {
    	dVdt_2 = 0.0;
    }

    // call the surrogate model for rheology
    // double shearRate = 0.3;

    // size_t temp_rheopos     = modena_model_inputs_argPos(rheologymodel, "temp");
    // size_t conv_rheopos     = modena_model_inputs_argPos(rheologymodel, "conv");
    // size_t shear_rheopos    = modena_model_inputs_argPos(rheologymodel, "shear");

    // modena_model_argPos_check(rheologymodel);

    // // set input vector
    // modena_inputs_set(inputs_rheo, temp_rheopos, T);
    // modena_inputs_set(inputs_rheo, conv_rheopos, EG_XOH);
    // modena_inputs_set(inputs_rheo, shear_rheopos, shearRate);

    // // call the model
    // int ret_rheo = modena_model_call (rheologymodel, inputs_rheo, outputs_rheo);

    // // terminate, if requested
    // if(ret_rheo != 0)
    // {
    //     modena_inputs_destroy (inputs_rheo);
    //     modena_outputs_destroy (outputs_rheo);
    //     modena_model_destroy (rheologymodel);
    // }

    // double mu_app = modena_outputs_get(outputs_rheo, 0);

    // Gelling point representation
    if(y[1] > 0.5)
    {
    	beta0	= 0.0;
    }
    else
    {
    	beta0	= beta0;
    }

    PDA(we, vi, mom, nNodes);
    growthSource(sgBA, sgCO2, we, vi, nNodes, mOrder, CO2_l, L_l, T, dVdt_2, dVdt_1);
    coalescenceSource(sc, we, vi, nNodes, mOrder, beta0);

    dydt[3] = -sgBA[1]*(p_1/(RR*y[2]))*(M_B/1000.0)*(1.0/rhoPolySurrgate);
	dydt[4]	=  sgBA[1]*(p_1/(RR*y[2]))*(M_B/1000.0)*(1.0/rhoPolySurrgate);
	dydt[6]	=  sgCO2[1]*(p_2/(RR*y[2]))*(M_CO2/1000.0)*(1.0/rhoPolySurrgate);
    dydt[5]	= -sgCO2[1]*(p_2/(RR*y[2]))*(M_CO2/1000.0)*(1.0/rhoPolySurrgate) + W_0*dydt[0]*(M_CO2/1000.0)*(1.0/rhoPolySurrgate);

	dydt[7] 	= sgBA[0] + sgCO2[0] + sc[0];
	dydt[8]  	= sgBA[1] + sgCO2[1] + sc[1];
	dydt[9]  	= sgBA[2] + sgCO2[2] + sc[2];
	dydt[10] 	= sgBA[3] + sgCO2[3] + sc[3];
}
/**
@fn main(int argc, char **argv)
@brief main function, initializes the state_type variables and performs the integration.

*/
int main(int argc, char **argv)
{
	readParams();
	// initial conditions
    state_type y(11);
    y[0]			= 0.0;
    y[1]			= 0.0;
    y[2]			= Temp0;
    y[3]			= L0;
    y[4]			= 1.0e-14; 
    y[5]			= 0.0;
    y[6]			= 1.0e-14;

    // moments initialization
    int nOfmom 		= 4;
    double momz[nOfmom];
    double pBA, pCO2, bubble_radius;
    mom_init(momz, init_size, nOfmom, sig, NN);

    y[7]			= momz[0];
    y[8]			= momz[1];
    y[9]			= momz[2];
    y[10]			= momz[3];
	double R = bubbleRadius(y[7], y[8]);
	air_g=y[8]/(1+y[8])*M_air*1e-3*(Pr+2*surfaceTension/R)/(RR*Temp0*rhoPoly);

	// initialize simpleKinetics variables
    // y[11]           = 5.0;
    // y[12]           = 5.0;
    // y[13]           = 0.2;
    // y[14]           = 0.0;
    // y[15]           = 0.0;
    // y[16]           = 0.0;
    // y[17]           = 0.0;
    // y[18]           = 0.0;
    // y[19]           = 300.0;


    runge_kutta4< state_type > stepper;

    for( double t=0.0 ; t<tend ; t+= dt )
    {
        /// @sa http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/steppers.html
		integrate_adaptive( make_controlled( abs_err , rel_err , error_stepper_type() ), QmomKinetics , y , t, t+dt , 1e-9 );
        write_kinetics(y, t);

        pBA           = partialPressureBA(y);
        pCO2          = partialPressureCO2(y);
        bubble_radius = bubbleRadius(y[7], y[8]);
        ddtpartialPressure(y, t, dt, dpdt, pOld, pBA, pCO2, bubble_radius);
    }
}

/* 

Different methods of integrations:

[ define_const_stepper
    runge_kutta4< state_type > stepper;
    integrate_const( stepper , harmonic_oscillator , x , 0.0 , 10.0 , 0.01 );
]

[ integrate_const_loop
    const double dt = 0.01;
    for( double t=0.0 ; t<10.0 ; t+= dt )
        stepper.do_step( harmonic_oscillator , x , t , dt );
 ]

 [ define_adapt_stepper
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
 ]

 [integrate_adapt_make_controlled
    integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-10 , 1.0e-6 ) ,
                        harmonic_oscillator , x , 0.0 , 10.0 , 0.01 );
 ]
[ integrate_const with abs and rel error
integrate_const( make_dense_output( 1.0e-6 , 1.0e-6 , runge_kutta_dopri5< state_type >() ) , sys , inout , t_start , t_end , dt );
irst two parameters are the absolute and the relative error tolerances
]
[ using adaptive integrate
double abs_err 	= 1.0e-12;
    double rel_err 	= 1.0e-10;
	integrate_adaptive( make_controlled< error_stepper_type >(abs_err , rel_err), kinetics, y, 0.0, 300.0, 0.01, write_kinetics );

]

*/
