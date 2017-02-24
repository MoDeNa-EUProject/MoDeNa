/*

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
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

Description
    Solving the two tank problem the fully integrated way.

Authors
    Henrik Rusche

Contributors
*/

#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

double flowRate
(
    const double D,
    const double rho0,
    const double p0,
    const double p2
)
{
    // from VDI Waermeatlas Lbd 4
    const double kappa = 1.4;

    const double etac = pow(2.0/(kappa+1.0), kappa/(kappa-1.0));

    double p1 = etac*p0;
    if(p2 > p1)
    {
        p1 = p2;
    }

    // for a nozzle
    const double Cd0 = 0.84;
    const double Cd1 = 0.66;
    // Not 100% sure this is correct - eqn. (10) in Lbd 5 uses misleading
    // nomenclature
    const double Cdg =
        Cd0 - Cd1*pow(p1/p0, 2.0) + (2*Cd1-Cd0)*pow(p1/p0, 3.0);

    const double Phi =
        sqrt
        (
            kappa/(kappa-1.0)
            *(pow(p1/p0, 2.0/kappa) - pow(p1/p0, (kappa+1.0)/kappa))
        );

    return M_PI*pow(D, 2.0)*Cdg*Phi*sqrt(2.0*rho0*p0);
}


int
main (int argc, char *argv[])
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

    while(t + deltat < tend + 1e-10)
    {
        t += deltat;

        if(p0 > p1)
        {
            double mdot = flowRate(D, rho0, p0, p1);
            m0 -= mdot*deltat;
            m1 += mdot*deltat;
        }
        else
        {
            double mdot = flowRate(D, rho0, p1, p0);
            m0 += mdot*deltat;
            m1 -= mdot*deltat;
        }

        rho0 = m0/V0;
        p0 = m0/V0*287.1*T;
        p1 = m1/V1*287.1*T;

        cout << "t = " << t << " p0 = " << p0 << " p1 = " << p1 << endl;
    }

    return 0;
}
