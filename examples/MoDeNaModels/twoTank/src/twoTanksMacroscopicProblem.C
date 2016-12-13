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
Solving the two tank problem the MoDeNa way.
A prototypical macros-scopic code embeds a micro-scale model (flowRate)
through the MoDeNa interface library.
@author     Henrik Rusche
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@ingroup    twoTank
*/

#include <stdio.h>
#include <iostream>
#include "modena.h"

using namespace std;

int
main(int argc, char *argv[])
{
    const double D = 0.01;

    double p0 = 3e5;
    double p1 = 10000;
    double V0 = 0.1;
    double V1 = 1;
    double T = 300;

    double t = 0.0;
    const double deltat = 1e-3;
    const double tend = 5.5;

    double m0 = p0*V0/287.1/T;
    double m1 = p1*V1/287.1/T;

    double rho0 = m0/V0;
    double rho1 = m1/V1;


    // Instantiate a model
    modena_model_t *model = modena_model_new("flowRate");

    if(modena_error_occurred())
    {
        return modena_error();
    }

    // Allocate memory and fetch arg positions
    modena_inputs_t *inputs = modena_inputs_new(model);
    modena_outputs_t *outputs = modena_outputs_new(model);

    size_t Dpos = modena_model_inputs_argPos(model, "D");
    size_t rho0Pos = modena_model_inputs_argPos(model, "rho0");
    size_t p0Pos = modena_model_inputs_argPos(model, "p0");
    size_t p1Byp0Pos = modena_model_inputs_argPos(model, "p1Byp0");

    modena_model_argPos_check(model);

    while(t + deltat < tend + 1e-10)
    {
        t += deltat;

        if(p0 > p1)
        {
            // Set input vector
            modena_inputs_set(inputs, Dpos, D);
            modena_inputs_set(inputs, rho0Pos, rho0);
            modena_inputs_set(inputs, p0Pos, p0);
            modena_inputs_set(inputs, p1Byp0Pos, p1/p0);

            // Call the model
            int ret = modena_model_call(model, inputs, outputs);

            // Terminate, if requested
            if(modena_error_occurred())
            {
                modena_inputs_destroy(inputs);
                modena_outputs_destroy(outputs);
                modena_model_destroy(model);

                return modena_error();
            }

            // Fetch result
            double mdot = modena_outputs_get(outputs, 0);

            m0 -= mdot*deltat;
            m1 += mdot*deltat;
        }
        else
        {
            // Set input vector
            modena_inputs_set(inputs, Dpos, D);
            modena_inputs_set(inputs, rho0Pos, rho1);
            modena_inputs_set(inputs, p0Pos, p1);
            modena_inputs_set(inputs, p1Byp0Pos, p0/p1);

            // Call the model
            int ret = modena_model_call(model, inputs, outputs);

            // Terminate, if requested
            if(modena_error_occurred())
            {
                modena_inputs_destroy(inputs);
                modena_outputs_destroy(outputs);
                modena_model_destroy(model);

                return modena_error();
            }

            // Fetch result
            double mdot = modena_outputs_get(outputs, 0);

            m0 += mdot*deltat;
            m1 -= mdot*deltat;
        }

        rho0 = m0/V0;
        rho1 = m1/V1;
        p0 = m0/V0*287.1*T;
        p1 = m1/V1*287.1*T;

        cout << "t = " << t << " rho0 = " << rho0 << " p0 = " << p0 << " p1 = " << p1 << endl;
    }

    modena_inputs_destroy(inputs);
    modena_outputs_destroy(outputs);
    modena_model_destroy(model);

    return 0;
}
