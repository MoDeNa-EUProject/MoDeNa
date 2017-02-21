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
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

Description
    This code calculates the flowRate through a nozzle as a function of the
    diameter, density and pressure upstream and pressure downstream. It uses
    expressions from VDI Waermeatlas Lbd 4

    In the simple twoTank example, this piece of code stands for a complex
    microscopic code - such as full 3D CFD simulation.

Authors
    Henrik Rusche

Contributors
*/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;


double flowRate
(
    const double D,
    const double rho0,
    const double p0,
    const double p2Byp0
)
{
    const double p2 = p0*p2Byp0;
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
    ifstream fi;
    fi.open ("in.txt");
    if(!fi.is_open())
    {
        cerr << "Could not open in.txt!" << endl;
        return 1;
    }

    double D, rho0, p0, p1Byp0;
    fi >> D >> rho0 >> p0 >> p1Byp0;

    double mdot = flowRate(D, rho0, p0, p1Byp0);

    double p1 = p0*p1Byp0;

    cout << "D = " << D
        << " rho0 = " << rho0
        << " p0 = " << p0
        << " p1/p0 = " << p1Byp0
        << " p1 = " << p1
        << " mdot = " << mdot
        << endl;

    ofstream fo;
    fo.open ("out.txt");

    fo << mdot << endl;

    fi.close();
    fo.close();
}

