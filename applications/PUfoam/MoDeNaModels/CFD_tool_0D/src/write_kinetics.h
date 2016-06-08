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
    out[0].open("./XW1.txt", std::ios::app);
    out[1].open("./XOH1.txt", std::ios::app);
    out[2].open("./T1.txt", std::ios::app);
    out[3].open("./L_l1.txt", std::ios::app);
    out[4].open("./L_g1.txt", std::ios::app);
    out[5].open("./CO2_l1.txt", std::ios::app);
    out[6].open("./CO2_g1.txt", std::ios::app);
    out[7].open("./m01.txt", std::ios::app);
    out[8].open("./m11.txt", std::ios::app);
    out[9].open("./m21.txt", std::ios::app);
    out[10].open("./m31.txt", std::ios::app);
    // out[11].open("./EG_XNCO1.txt", std::ios::app);
    // out[12].open("./EG_XOH1.txt", std::ios::app);
    // out[13].open("./XH2O1.txt", std::ios::app);
    // out[14].open("./CO21.txt", std::ios::app);
    // out[15].open("./PENTANE1.txt", std::ios::app);
    // out[16].open("./POLYMER1.txt", std::ios::app);
    // out[17].open("./POLYMERBLOW1.txt", std::ios::app);
    // out[18].open("./UREA1.txt", std::ios::app);
    // out[19].open("./R_1_temp1.txt", std::ios::app);

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
            // Calling the model for density reaction mixture
            size_t T_denpos     = modena_model_inputs_argPos(density_reaction_mixturemodel, "T");
            size_t XOH_denpos   = modena_model_inputs_argPos(density_reaction_mixturemodel, "XOH");
            modena_model_argPos_check(density_reaction_mixturemodel);
            // set input vector
            modena_inputs_set(inputs_den, T_denpos, y[2]);
            modena_inputs_set(inputs_den, XOH_denpos, y[1]);
            // call the model
            int ret_den = modena_model_call (density_reaction_mixturemodel, inputs_den, outputs_den);
            // terminate, if requested
            if(ret_den != 0)
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
            modena_inputs_set(inputs_den, T_denpos, y[2]);
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
        default:
            cerr << "Invalid density model" << endl;
            exit(1);
    }

    double p1,p2;
    p1 = partialPressureBA(y);
    p2 = partialPressureCO2(y);
	double rho_bubble 	= ((p1+p2)/(RR*y[2]))*(y[6]*M_CO2 + y[4]*M_B)/(fmax((1000.0*(y[6] + y[4])),1.0e-8));
	// double rho_foam 	= (rho_bubble*(y[8]/(1.0+y[8])) + (1.0+L0)*rhoPolySurrgate*(1.0 - (y[8]/(1.0+y[8]))));
    double rho_foam 	= (rho_bubble*(y[8]/(1.0+y[8])) + rhoPolySurrgate*(1.0 - (y[8]/(1.0+y[8]))));

    // print out the strut content
    // set input vector
    modena_inputs_set(inputs_strutContent, rho_foam_Pos, rho_foam);

    // call the model
    int ret_strutContent = modena_model_call (strutContentmodel, inputs_strutContent, outputs_strutContent);
    if ((tend - t) < 2)
    {
        cout << "final foam density: " << rho_foam << endl;
        cout << "strut content: " << modena_outputs_get(outputs_strutContent, 0) << endl;
    }

    double thermalConductivity;
    if (rho_foam > 48.0)
    {
        thermalConductivity = 8.7006e-8*rho_foam*rho_foam + 8.4674e-5*rho_foam
                             + 1.16e-2;
    }
    else
    {
        thermalConductivity = 9.3738e-6*rho_foam*rho_foam - 7.3511e-4*rho_foam
                             + 2.956e-2;
    }
    // surrogate model for thermal conductivity
    modena_inputs_set(inputs_thermalConductivity, porosity_Pos, (1.0 - rho_foam/rhoPolySurrgate));
    double R = bubbleRadius(y[7], y[8]);
    modena_inputs_set(inputs_thermalConductivity, cell_size_Pos, (2.0*R));
    modena_inputs_set(inputs_thermalConductivity, temp_Pos, y[2]);
    modena_inputs_set(inputs_thermalConductivity, X_CO2_Pos, (p2/(p1+p2)));
    modena_inputs_set(inputs_thermalConductivity, X_O2_Pos, 0.0);
    modena_inputs_set(inputs_thermalConductivity, X_N2_Pos, 0.0);
    modena_inputs_set(inputs_thermalConductivity, X_Cyp_Pos, (p1/(p1+p2)));
    double st_c;
    st_c = modena_outputs_get(outputs_strutContent, 0);
    modena_inputs_set(inputs_thermalConductivity, strut_c_Pos, st_c);
    int ret_thermalConductivitymodel = modena_model_call (thermalConductivitymodel, inputs_thermalConductivity, outputs_thermalConductivity);
    if(modena_error_occurred())
    {
        exit(modena_error());
    }
    double the_con;
    the_con = modena_outputs_get(outputs_thermalConductivity, 0);
    ofstream thermalOut;
    thermalOut.open("./thermalConductivity.txt", std::ios::app);
    thermalOut << the_con << '\t' << thermalConductivity << '\t' << rho_foam << '\n';
    thermalOut.close();


	ofstream rho_bubbleout;
	rho_bubbleout.open("./rho_bubble.txt", std::ios::app);
	rho_bubbleout << t << '\t' << rho_bubble << '\n';
	rho_bubbleout.close();

	ofstream rho_foamout;
	rho_foamout.open("./rho_foam.txt", std::ios::app);
	rho_foamout << t << '\t' << rho_foam << '\n';
	rho_foamout.close();

    ofstream R_bubble;
    R_bubble.open("./R_bubble.txt", std::ios::app);
    R_bubble << t << '\t' << R << '\n';
    R_bubble.close();

    momentsConverter(y,t);

	double p_1 = partialPressureBA(y);
	ofstream p1out;
	p1out.open("./p_1.txt", std::ios::app);
	p1out << t << '\t' << p_1 << '\n';
	p1out.close();

	double p_2 = partialPressureCO2(y);
	ofstream p2out;
	p2out.open("./p_2.txt", std::ios::app);
	p2out << t << '\t' << p_2 << '\n';
	p2out.close();
}
