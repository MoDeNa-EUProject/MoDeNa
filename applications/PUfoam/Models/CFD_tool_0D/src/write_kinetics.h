/** @file write_kinetics.h
    @brief write the results into text files
    @fn void write_kinetics( const state_type &y , const double t )
    @param const state_type &y -  vector of all the variables
    @param const double t - time
    @return void
*/
void write_kinetics( const state_type &y , const double t );

void write_kinetics( const state_type &y , const double t )
{
// write_kinetics - write the results into text files
// @param - const state_type &y -  vector of all the variables
// @param - const double t - time
    int n;
    n = y.size();

    ofstream out[n];
    out[0].open("../results/XW1.txt", std::ios::app);
    out[1].open("../results/XOH1.txt", std::ios::app);
    out[2].open("../results/T1.txt", std::ios::app);
    out[3].open("../results/L_l1.txt", std::ios::app);
    out[4].open("../results/L_g1.txt", std::ios::app);
    out[5].open("../results/CO2_l1.txt", std::ios::app);
    out[6].open("../results/CO2_g1.txt", std::ios::app);
    out[7].open("../results/m01.txt", std::ios::app);
    out[8].open("../results/m11.txt", std::ios::app);
    out[9].open("../results/m21.txt", std::ios::app);
    out[10].open("../results/m31.txt", std::ios::app);
    // out[11].open("../results/EG_XNCO1.txt", std::ios::app);
    // out[12].open("../results/EG_XOH1.txt", std::ios::app);
    // out[13].open("../results/XH2O1.txt", std::ios::app);
    // out[14].open("../results/CO21.txt", std::ios::app);
    // out[15].open("../results/PENTANE1.txt", std::ios::app);
    // out[16].open("../results/POLYMER1.txt", std::ios::app);
    // out[17].open("../results/POLYMERBLOW1.txt", std::ios::app);
    // out[18].open("../results/UREA1.txt", std::ios::app);
    // out[19].open("../results/R_1_temp1.txt", std::ios::app);

    for (int i = 0; i < n; i++)
    {
        if(i == 11)
        {
            out[i] << t << '\t' << (1.0-y[i]/5.0) << '\n';
        }
        else if (i == 12)
        {
            out[i] << t << '\t' << (1.0-y[i]/5.0) << '\n';
        }
        else if (i == 13)
        {
            out[i] << t << '\t' << (1.0-y[i]/0.2) << '\n';
        }
        else
        {
            out[i] << t << '\t' << y[i] << '\n';
        }
        out[i].close();
    }

    double rhoPolySurrgate;
    switch (denMod) {
		case 1: {
        	// // Calling the model for density reaction mixture
         //    size_t T_denpos     = modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
         //    size_t XOH_denpos   = modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
         //    modena_model_argPos_check(density_reaction_mixturemodel);
         //    // set input vector
         //    double EG_XOH;
         //    double init_EG_OH = 5.0;
         //    EG_XOH      = 1.0 - (y[12]/init_EG_OH);
         //    modena_inputs_set(inputs_den, T_denpos, y[2]);
         //    modena_inputs_set(inputs_den, XOH_denpos, EG_XOH);
         //    // call the model
         //    int ret_den = modena_model_call (density_reaction_mixturemodel, inputs_den, outputs_den);
         //    // terminate, if requested
         //    if(ret_den != 0)
         //    {
         //        modena_inputs_destroy (inputs_den);
         //        modena_outputs_destroy (outputs_den);
         //        modena_model_destroy (density_reaction_mixturemodel);
         //        exit(ret_den);
         //        //return ret_den;
         //    }
         //    rhoPolySurrgate = modena_outputs_get(outputs_den, 0);
        }
        case 2:
            rhoPolySurrgate = rhoPoly;
    }

    double p1,p2;
    p1 = partialPressureBA(y);
    p2 = partialPressureCO2(y);
	double rho_bubble 	= ((p1+p2)/(RR*y[2]))*(y[6]*M_CO2 + y[4]*M_B)/(fmax((1000.0*(y[6] + y[4])),1.0e-8));
	// double rho_foam 	= (rho_bubble*(y[8]/(1.0+y[8])) + (1.0+L0)*rhoPolySurrgate*(1.0 - (y[8]/(1.0+y[8]))));
    double rho_foam 	= (rho_bubble*(y[8]/(1.0+y[8])) + rhoPolySurrgate*(1.0 - (y[8]/(1.0+y[8]))));

	ofstream rho_bubbleout;
	rho_bubbleout.open("../results/rho_bubble.txt", std::ios::app);
	rho_bubbleout << t << '\t' << rho_bubble << '\n';
	rho_bubbleout.close();

	ofstream rho_foamout;
	rho_foamout.open("../results/rho_foam.txt", std::ios::app);
	rho_foamout << t << '\t' << rho_foam << '\n';
	rho_foamout.close();

    ofstream R_bubble;
    double R = bubbleRadius(y[7], y[8]);
    R_bubble.open("../results/R_bubble.txt", std::ios::app);
    R_bubble << t << '\t' << R << '\n';
    R_bubble.close();

    momentsConverter(y,t);

	double p_1 = partialPressureBA(y);
	ofstream p1out;
	p1out.open("../results/p_1.txt", std::ios::app);
	p1out << t << '\t' << p_1 << '\n';
	p1out.close();

	double p_2 = partialPressureCO2(y);
	ofstream p2out;
	p2out.open("../results/p_2.txt", std::ios::app);
	p2out << t << '\t' << p_2 << '\n';
	p2out.close();
}
