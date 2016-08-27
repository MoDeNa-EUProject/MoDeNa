/** @cond

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2016 MoDeNa Consortium, All rights reserved.

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
@brief      macro-scale tool for the foaming process.
@authors    Mohsen Karimi, Daniele Marchisio, Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@ingroup    app_foaming
*/

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/range/numeric.hpp>
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
// #include "bubbleRadius.h"
#include "partialPressure.h"
/**
@var dpdt[2]: global double array
@brief This is used to compute the partial pressures.

@var pOld[2]: global double array variable
@brief This is to hold the old pressure values during the partial pressure calculations.
*/
double dpdt[2] = {};
double pOld[2] = {};

#include "modenaData.h"
#include "momentsConverter.h"
#include "write_kinetics.h"
#include "determinant.h"
#include "HankelHadamard.h"
#include "differenceTable.h"
#include "McGrawCorrection.h"
#include "WrightCorrection.h"
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
    // ---RF-1-private---
    // dydt[11]: Catalyst_1
    // dydt[12]: CE_A0
    // dydt[13]: CE_A1
    // dydt[14]: CE_B
    // dydt[15]: CE_B2
    // dydt[16]: CE_I0
    // dydt[17]: CE_I1
    // dydt[18]: CE_I2
    // dydt[19]: CE_PBA
    // dydt[20]: CE_Breac
    // dydt[21]: CE_Areac0
    // dydt[22]: CE_Areac1
    // dydt[23]: CE_Ireac0
    // dydt[24]: CE_Ireac1
    // dydt[25]: CE_Ireac2
    // dydt[26]: Bulk
    // dydt[27]: R_1
    // dydt[28]: R_1_mass
    // dydt[29]: R_1_temp
    // dydt[30]: R_1_vol

	int nNodes = 2;
	int mOrder[6] = {0, 1, 2, 3, 4, 5};
	double mom[2*nNodes], we[nNodes], vi[nNodes], sgBA[2*nNodes], sgCO2[2*nNodes], sc[2*nNodes];
	double L_l, L_g, CO2_l, CO2_g, T, Lm, c1;
	double rhoPolySurrgate;
	double beta0 	= 0.0;
    double XW, XOH;

    // RF-1 variables
    double Catalyst_1, CE_A0, CE_A1, CE_B, CE_B2,
           CE_I0, CE_I1, CE_I2, CE_PBA, CE_Breac,
           CE_Areac0, CE_Areac1, CE_Ireac0, CE_Ireac1,
           CE_Ireac2, Bulk, R_1, R_1_mass, R_1_temp, R_1_vol;

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
    Catalyst_1  = y[11];
    CE_A0       = y[12];
    CE_A1       = y[13];
    CE_B        = y[14];
    CE_B2       = y[15];
    CE_I0       = y[16];
    CE_I1       = y[17];
    CE_I2       = y[18];
    CE_PBA      = y[19];
    CE_Breac    = y[20];
    CE_Areac0   = y[21];
    CE_Areac1   = y[22];
    CE_Ireac0   = y[23];
    CE_Ireac1   = y[24];
    CE_Ireac2   = y[25];
    Bulk        = y[26];
    R_1         = y[27];
    R_1_mass    = y[28];
    R_1_temp    = y[29];
    R_1_vol     = y[30];

    double Rx;
    // Calling the RF-1 model
    switch (kinMod)
    {
        case 1:
            Rx=1;
            break;
        case 2:
            if (y[1]<X_gel) {
                Rx=1;
            } else if (y[1]<0.87) {
                Rx=-2.027*y[1]+2.013;
            } else {
                Rx=3.461*y[1]-2.761;
            }
            break;
        case 3:
            // set input vector
            modena_inputs_set(inputs_kinetics, kineticTime_Pos, t);
            modena_inputs_set(inputs_kinetics, Catalyst_1_Pos, Catalyst_1);
            modena_inputs_set(inputs_kinetics, CE_A0_Pos, CE_A0);
            modena_inputs_set(inputs_kinetics, CE_A1_Pos, CE_A1);
            modena_inputs_set(inputs_kinetics, CE_B_Pos, CE_B);
            modena_inputs_set(inputs_kinetics, CE_B2_Pos, CE_B2);
            modena_inputs_set(inputs_kinetics, CE_I0_Pos, CE_I0);
            modena_inputs_set(inputs_kinetics, CE_I1_Pos, CE_I1);
            modena_inputs_set(inputs_kinetics, CE_I2_Pos, CE_I2);
            modena_inputs_set(inputs_kinetics, CE_PBA_Pos, CE_PBA);
            modena_inputs_set(inputs_kinetics, CE_Breac_Pos, CE_Breac);
            modena_inputs_set(inputs_kinetics, CE_Areac0_Pos, CE_Areac0);
            modena_inputs_set(inputs_kinetics, CE_Areac1_Pos, CE_Areac1);
            modena_inputs_set(inputs_kinetics, CE_Ireac0_Pos, CE_Ireac0);
            modena_inputs_set(inputs_kinetics, CE_Ireac1_Pos, CE_Ireac1);
            modena_inputs_set(inputs_kinetics, CE_Ireac2_Pos, CE_Ireac2);
            modena_inputs_set(inputs_kinetics, Bulk_Pos, Bulk);
            modena_inputs_set(inputs_kinetics, R_1_Pos, R_1);
            modena_inputs_set(inputs_kinetics, R_1_mass_Pos, R_1_mass);
            modena_inputs_set(inputs_kinetics, R_1_temp_Pos, (T-273.15));
            modena_inputs_set(inputs_kinetics, R_1_vol_Pos, R_1_vol);

            // call the model
            int ret_kinetics = modena_model_call(kinetics, inputs_kinetics, outputs_kinetics);

            // terminate, if requested
            if(modena_error_occurred())
            {
                modena_inputs_destroy (inputs_kinetics);
                modena_outputs_destroy (outputs_kinetics);
                modena_model_destroy (kinetics);
                cout << "Modena Error:" << (modena_error()) << endl;
            }

         // get the source terms for simpleKinetics
            dydt[11] = modena_outputs_get(outputs_kinetics, source_Catalyst_1_Pos);
            dydt[12] = modena_outputs_get(outputs_kinetics, source_CE_A0_Pos);
            dydt[13] = modena_outputs_get(outputs_kinetics, source_CE_A1_Pos);
            dydt[14] = modena_outputs_get(outputs_kinetics, source_CE_B_Pos);
            dydt[15] = modena_outputs_get(outputs_kinetics, source_CE_B2_Pos);
            dydt[16] = modena_outputs_get(outputs_kinetics, source_CE_I0_Pos);
            dydt[17] = modena_outputs_get(outputs_kinetics, source_CE_I1_Pos);
            dydt[18] = modena_outputs_get(outputs_kinetics, source_CE_I2_Pos);
            dydt[19] = modena_outputs_get(outputs_kinetics, source_CE_PBA_Pos);
            dydt[20] = modena_outputs_get(outputs_kinetics, source_CE_Breac_Pos);
            dydt[21] = modena_outputs_get(outputs_kinetics, source_CE_Areac0_Pos);
            dydt[22] = modena_outputs_get(outputs_kinetics, source_CE_Areac1_Pos);
            dydt[23] = modena_outputs_get(outputs_kinetics, source_CE_Ireac0_Pos);
            dydt[24] = modena_outputs_get(outputs_kinetics, source_CE_Ireac1_Pos);
            dydt[25] = modena_outputs_get(outputs_kinetics, source_CE_Ireac2_Pos);
            dydt[26] = modena_outputs_get(outputs_kinetics, source_Bulk_Pos);
            dydt[27] = modena_outputs_get(outputs_kinetics, source_R_1_Pos);
            dydt[28] = modena_outputs_get(outputs_kinetics, source_R_1_mass_Pos);
            dydt[29] = modena_outputs_get(outputs_kinetics, source_R_1_temp_Pos);
            dydt[30] = modena_outputs_get(outputs_kinetics, source_R_1_vol_Pos);

			if (W_0 > 1e-8) {
				dydt[0] = -dydt[15]/(W_0/1000);
			} else {
				dydt[0] = 0;
			}
            dydt[1] = -(dydt[12] + dydt[13])/(OH_0)*1000;
            break;
    }

    switch (denMod)
    {
        case 1:
        {
            // Calling the model for density reaction mixture
            size_t T_denpos     = modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
            size_t XOH_denpos   = modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
            modena_model_argPos_check(density_reaction_mixturemodel);

            // // set input vector
            modena_inputs_set(inputs_den, T_denpos, T);
            modena_inputs_set(inputs_den, XOH_denpos, XOH);

            // // call the model
            int ret_den = modena_model_call (density_reaction_mixturemodel, inputs_den, outputs_den);

            if (ret_den != 0)
            {
                modena_inputs_destroy (inputs_den);
                modena_outputs_destroy (outputs_den);
                modena_model_destroy (density_reaction_mixturemodel);
                exit(ret_den);
            }

            rhoPolySurrgate = modena_outputs_get(outputs_den, 0);
            break;
        }
        case 2:
            rhoPolySurrgate = rhoPoly;
            break;
		case 3:
		{
			// Calling the PCSAFT model for density reaction mixture
            size_t T_denpos     = modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
            modena_model_argPos_check(density_reaction_mixturemodel);
            modena_inputs_set(inputs_den, T_denpos, T);
            // // call the model
            int ret_den = modena_model_call (density_reaction_mixturemodel, inputs_den, outputs_den);
            if (ret_den != 0)
            {
                modena_inputs_destroy (inputs_den);
                modena_outputs_destroy (outputs_den);
                modena_model_destroy (density_reaction_mixturemodel);
                exit(ret_den);
            }
            rhoPolySurrgate = modena_outputs_get(outputs_den, 0);
            break;
		}
    }
    if (kinMod == 1 || kinMod == 2)
    {
        // ODEs
        dydt[0]     = A_W*exp(-E_W/(RR*y[2]))*(1-y[0]);
        if(dydt[0] < 0.0)
        {
            dydt[0] = 0.0;
        }
        if(W_0 < 0.0)
        {
            dydt[0] = 0.0;
        }

        dydt[1]     = Rx*A_OH*exp(-E_OH/(RR*y[2]))*OH_0*(1-y[1])*(NCO_0/OH_0 - 2.0*y[0]*W_0/OH_0 - y[1]);
        if(dydt[1] < 0.0)
        {
            dydt[1] = 0.0;
        }
        if (dilution) {
            dydt[0]=dydt[0]/(1+y[3]*rhoPolySurrgate/rhoBL);
            dydt[1]=dydt[1]/(1+y[3]*rhoPolySurrgate/rhoBL);
        }
    }

    // call the surrogate model for rheology
    if (apparentViscosity)
    {
        double shearRate = 0.05;
        // // set input vector
        modena_inputs_set(inputs_rheo, temp_rheopos, T);
        modena_inputs_set(inputs_rheo, conv_rheopos, XOH);
        modena_inputs_set(inputs_rheo, shear_rheopos, shearRate);
		modena_inputs_set(inputs_rheo, m0_rheopos, mom[0]);
		modena_inputs_set(inputs_rheo, m1_rheopos, mom[1]);
        // // call the model
        int ret_rheo = modena_model_call (rheologymodel, inputs_rheo, outputs_rheo);
        // // terminate, if requested
        if(modena_error_occurred())
        {
            modena_inputs_destroy (inputs_rheo);
            modena_outputs_destroy (outputs_rheo);
            modena_model_destroy (rheologymodel);
            cout << "Modena Error: " << (modena_error()) << endl;
        }
        double mu_app = modena_outputs_get(outputs_rheo, 0);
        // cout << "apparent viscosity: " << mu_app << endl;
    }

    // Gelling point representation
    if(y[1] > X_gel)
    {
        beta0   = 0.0;
    }
    else
    {
        beta0   = beta0;
    }

    if (realizabilityCheck)
    {
        cout << "----Step 1----" << endl;
        cout << "moments before check realizability: " << endl;
        printMoms(mom,nNodes);

        // Check realizability
        cout << "----Step 2----" << endl;
        cout << "check positivity" << endl;
        if (!momentsPositivite(mom,nNodes))
        {
            double initial_mom[2*nNodes];
            int nMoms = 2*nNodes;
            mom_init(initial_mom, init_size, nMoms, sig, NN);

        for (int i = 0; i < 2*nNodes; i++)
        {
            if (mom[i] < 0.0)
            {
               mom[i] = initial_mom[i];
            }
        }
        }
        printMoms(mom,nNodes);


        int realizable = 0;
        cout << "----Step 3----" << endl;
        cout << "first Hankel-Hadamard check" << endl;
        realizable = HankelHadamard(mom, nNodes);
        cout << "Hankel-Hadamard check (0:realizable, 1:unrealizable) --> " << realizable <<  endl;

        if (realizable == 1)
        {
            // double M0 = mom[0];
            // normalizeMom(mom, nNodes);
            cout << "----Step 4----" << endl;
            cout << "McGrawCorrection" << endl;
            McGrawCorrection(mom, nNodes);
            cout << "moments after McGrawCorrection: " << endl;
            printMoms(mom,nNodes);
            // denormalizeMom(mom, M0, nNodes);
        }
        cout << "----Step 5----" << endl;
        cout << "second Hankel-Hadamard check" << endl;
        realizable = HankelHadamard(mom, nNodes);
        cout << "Hankel-Hadamard check (0:realizable, 1:unrealizable) --> " << realizable <<  endl;

        if (realizable == 1)
        {
            cout << "----Step 6----" << endl;
            cout << "WrightCorrection" << endl;
            WrightCorrection(mom, nNodes);
            cout << "moments after WrightCorrection: " << endl;
            printMoms(mom,nNodes);
        }
        // cout << "moments after check realizability: " << endl;
        // printMoms(mom,nNodes);
    }

    PDA(we, vi, mom, nNodes);

	Lm 			= LMax(T);
    // calling the surogate models for bubble growth rates.

    // partial pressure within bubbles due to the evaporation of physical blowing agent
    double p_1  = partialPressureBA(y);
    // // partial pressure within bubbles due to the generation of CO2
    double p_2  = partialPressureCO2(y);
    double c_1  = L_l*rhoPolySurrgate*1000.0/M_B;
    double c_2  = CO2_l*rhoPolySurrgate*1000.0/M_CO2;

    if (bubbleMode == "two nodes")
    {
        double radiusGrowthBA[nNodes],
            radiusGrowthCO2[nNodes],
            volumeGrowthBA[nNodes],
            volumeGrowthCO2[nNodes],
            nodeRadii[nNodes];
        int     ret_bblgr1[nNodes], ret_bblgr2[nNodes];
        for (int i = 0; i < nNodes; i++)
        {
            nodeRadii[i] = nodeRadius(vi[i]);
            // set input vector
            modena_inputs_set(inputs_bblgr1, Tbblgr1pos, T);
            modena_inputs_set(inputs_bblgr1, Rbblgr1pos, nodeRadii[i]);
            // modena_inputs_set(inputs_bblgr1, KH1bblgr1pos, KH1);
            modena_inputs_set(inputs_bblgr1, c_1bblgr1pos, c_1);
            modena_inputs_set(inputs_bblgr1, p_1bblgr1pos, p_1);
            // call the bblgr1 model
            ret_bblgr1[i]       = modena_model_call (bblgr1, inputs_bblgr1, outputs_bblgr1);
            radiusGrowthBA[i]   = modena_outputs_get(outputs_bblgr1, 0);
            // set input vector
            modena_inputs_set(inputs_bblgr2, Tbblgr2pos, T);
            modena_inputs_set(inputs_bblgr2, Rbblgr2pos, nodeRadii[i]);
            // modena_inputs_set(inputs_bblgr2, KH2bblgr2pos, KH2);
            modena_inputs_set(inputs_bblgr2, c_2bblgr2pos, c_2);
            modena_inputs_set(inputs_bblgr2, p_2bblgr2pos, p_2);
            // call the bblgr2 model
            ret_bblgr2[i]       = modena_model_call (bblgr2, inputs_bblgr2, outputs_bblgr2);
            // cout << "ret_bblgr2[" << i << "] = " << ret_bblgr2[i] << endl;
            radiusGrowthCO2[i]  = modena_outputs_get(outputs_bblgr2, 0);
            // cout << "radiusGrowthCO2[" << i << "] = " << radiusGrowthCO2[i] << endl;

            if(modena_error_occurred())
            {
                cout << modena_error() << endl;
            }
            volumeGrowthBA[i]   = (radiusGrowthBA[i]*RR*T)/(p_1);
            volumeGrowthCO2[i]  = (radiusGrowthCO2[i]*RR*T)/(p_2);

            if (volumeGrowthBA[i] < 0.0 || radiusGrowthBA[i] < 0.0 || L0 < 1.0e-8 || y[1] > X_gel)
            {
                volumeGrowthBA[i]   = 0.0;
            }
            if (volumeGrowthCO2[i] < 0.0 || radiusGrowthCO2[i] < 0.0 || W_0 < 1.0e-8 || y[1] > X_gel)
            {
                volumeGrowthCO2[i]  = 0.0;
            }
        }

        growthSource(sgBA, sgCO2, we, vi, nNodes, mOrder, CO2_l, L_l, T, volumeGrowthBA, volumeGrowthCO2);
    }
    else if (bubbleMode == "mean radius")
    {
        // bubble radius for bblgr1 and bblgr2 model
        double R    = bubbleRadius(mom[0], mom[1]);
        // set input vector
        modena_inputs_set(inputs_bblgr1, Tbblgr1pos, T);
        modena_inputs_set(inputs_bblgr1, Rbblgr1pos, R);
        // modena_inputs_set(inputs_bblgr1, KH1bblgr1pos, KH1);
        modena_inputs_set(inputs_bblgr1, c_1bblgr1pos, c_1);
        modena_inputs_set(inputs_bblgr1, p_1bblgr1pos, p_1);
        // set input vector
        modena_inputs_set(inputs_bblgr2, Tbblgr2pos, T);
        modena_inputs_set(inputs_bblgr2, Rbblgr2pos, R);
        // modena_inputs_set(inputs_bblgr2, KH2bblgr2pos, KH2);
        modena_inputs_set(inputs_bblgr2, c_2bblgr2pos, c_2);
        modena_inputs_set(inputs_bblgr2, p_2bblgr2pos, p_2);

        // call the bblgr1 model
        int ret_bblgr_1 = modena_model_call (bblgr1, inputs_bblgr1, outputs_bblgr1);
        // call the bblgr2 model
        int ret_bblgr_2 = modena_model_call (bblgr2, inputs_bblgr2, outputs_bblgr2);

        double G1, G2;
        double dVdt_1[nNodes], dVdt_2[nNodes];
        for (int i = 0; i < nNodes; i++)
        {
            dVdt_1[i] = 0.0;
            dVdt_2[i] = 0.0;
        }
        G1 = modena_outputs_get(outputs_bblgr1, 0);
        G2 = modena_outputs_get(outputs_bblgr2, 0);
        dVdt_1[0] = (G1*RR*T)/(p_1);
        if (dVdt_1[0] < 0.0 || G1 < 0.0 || L0<1e-8 || y[1]>X_gel)
        {
            dVdt_1[0] = 0.0;
        }
        dVdt_1[1] = dVdt_1[0];
        dVdt_2[0] = (G2*RR*T)/(p_2);
        if (dVdt_2[0] < 0.0 || G2 < 0.0 || W_0<1e-8 || y[1]>X_gel)
        {
            dVdt_2[0] = 0.0;
        }
        dVdt_2[1] = dVdt_2[0];

        growthSource(sgBA, sgCO2, we, vi, nNodes, mOrder, CO2_l, L_l, T, dVdt_1, dVdt_2);
    }
    else
    {
        cerr << "Invalid choice of bubbleMode!" << endl;
        exit(1);
    }

    coalescenceSource(sc, we, vi, nNodes, mOrder, beta0);

    dydt[3] = -sgBA[1]*(p_1/(RR*y[2]))*(M_B/1000.0)*(1.0/rhoPolySurrgate);
	dydt[4]	=  sgBA[1]*(p_1/(RR*y[2]))*(M_B/1000.0)*(1.0/rhoPolySurrgate);
	dydt[6]	=  sgCO2[1]*(p_2/(RR*y[2]))*(M_CO2/1000.0)*(1.0/rhoPolySurrgate);
    dydt[5]	= -sgCO2[1]*(p_2/(RR*y[2]))*(M_CO2/1000.0)*(1.0/rhoPolySurrgate) +\
	 			W_0*dydt[0]*(M_CO2/1000.0)*(1.0/rhoPolySurrgate);

	dydt[7] 	= sgBA[0] + sgCO2[0] + sc[0];
	dydt[8]  	= sgBA[1] + sgCO2[1] + sc[1];
	dydt[9]  	= sgBA[2] + sgCO2[2] + sc[2];
	dydt[10] 	= sgBA[3] + sgCO2[3] + sc[3];
	// temperature
    C_TOT 		= C_Poly + CO2_g*C_CO2 + L_g*C_BG + L_l*C_BL;
    dydt[2] 	= (-DH_OH*OH_0)/(rhoPolySurrgate*C_TOT)*dydt[1]+\
				  (-DH_W*W_0)/(rhoPolySurrgate*C_TOT)*dydt[0]+\
				  lambda/C_TOT*dydt[3];
}
/**
@fn main(int argc, char **argv)
@brief main function, initializes the state_type variables and performs the integration.

*/
int main(int argc, char **argv)
{
	readParams();
    #include "modenaCalls.h"

	// initial conditions
    state_type y(31);
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

    if(momentsPositivite(momz,nOfmom))
    {
        y[7]            = momz[0];
        y[8]            = momz[1];
        y[9]            = momz[2];
        y[10]           = momz[3];
    }
	double R = bubbleRadius(y[7], y[8]);
	air_g=y[8]/(1+y[8])*M_air*1e-3*(Pr+2*surfaceTension/R)/(RR*Temp0*rhoPoly);

    // initialize RF-1 variables
    y[11]           = catalyst*1e-3;
    y[12]           = polyol1_ini*1e-3;
    y[13]           = polyol2_ini*1e-3;
    y[14]           = amine_ini*1e-3;
	y[15]           = max(W_0*1e-3,0.0);
    y[16]           = isocyanate1_ini*1e-3;
    y[17]           = isocyanate2_ini*1e-3;
    y[18]           = isocyanate3_ini*1e-3;
    y[19]           = 0.0;
    y[20]           = 0.0;
    y[21]           = 0.0;
    y[22]           = 0.0;
    y[23]           = 0.0;
    y[24]           = 0.0;
    y[25]           = 0.0;
    y[26]           = 0.0;
    y[27]           = 0.0;
    y[28]           = 0.0;
    y[29]           = 0.0;
    y[30]           = 0.0;


    runge_kutta4< state_type > stepper;

    ofstream file;
    file.open("resultsKinetics.txt");
    file.setf(ios::scientific | ios::showpoint);
    cout.precision(20);
    cout.setf(ios::fixed | ios::showpoint);

    file << setw(12) << "t" << setw(12) << "Catalyst_1"
                            << setw(12) << "CE_A0"
                            << setw(12) << "CE_A1"
                            << setw(12) << "CE_B"
                            << setw(12) << "CE_B2"
                            << setw(12) << "CE_I0"
                            << setw(12) << "CE_I1"
                            << setw(12) << "CE_I2"
                            << setw(12) << "CE_PBA"
                            << setw(12) << "CE_Breac"
                            << setw(12) << "CE_Areac0"
                            << setw(12) << "CE_Areac1"
                            << setw(12) << "CE_Ireac0"
                            << setw(12) << "CE_Ireac1"
                            << setw(12) << "CE_Ireac2"
                            << endl;

    for( double t=0.0 ; t<tend ; t+= dt )
    {
        cout << "Time = " << t << endl;
        /// @sa http://headmyshoulder.github.io/odeint-v2/doc/boost_numeric_odeint/odeint_in_detail/steppers.html
		integrate_adaptive( make_controlled( abs_err , rel_err , error_stepper_type() ), QmomKinetics , y , t, t+dt , 1e-9 );
        file << setw(12) << t << " " << setw(12) << y[11] << " " << setw(12) << y[12] << " " << setw(12) << y[13] << " "
             << setw(12) << y[14] << " " << setw(12) << y[15] << " " << setw(12) << y[16] << " "
             << setw(12) << y[17] << " " << setw(12) << y[18] << " " << setw(12) << y[19] << " "
             << setw(12) << y[20] << " " << setw(12) << y[21] << " " << setw(12) << y[22] << " "
             << setw(12) << y[23] << " " << setw(12) << y[24] << " " << setw(12) << y[25] << " "
             << endl;
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
