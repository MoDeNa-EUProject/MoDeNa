'''@cond

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
@endcond'''

"""
@file
Surrogate function and model definitions for thermal conductivity of mixture of
blowing agents.

@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
from modena import *
from fireworks.utilities.fw_utilities import explicit_serialize
import gasConductivity

## simple weighted average
weightedAverageCode=r'''
#include "modena.h"
#include "math.h"
#include "stdio.h"

void gasMixtureConductivity
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    double kgasmix=0; // gas mixture conductivity
    double sumx=0;
    int i;
    for (i=0;i<x_size;i++) {
        sumx = sumx + x[i];
    }
    for (i=0;i<x_size;i++) {
        kgasmix = kgasmix + gas_thermal_conductivity[i]*x[i]/sumx;
    }
    outputs[0] = kgasmix;
}
'''
## Dohrn model
#
#  see [link](http://dx.doi.org/10.1016/j.fluid.2007.07.059)
dohrnCode=r'''
#include "modena.h"
#include "math.h"
#include "stdio.h"

void gasMixtureConductivity
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    double kgasmix=0; // gas mixture conductivity
    double Rg=8.314; // gas constant
    double y=0;
    double xi[4]; // normalized molar fractions
    double k[4]; // thermal conductivities
    double Tc[4]={304.17,511.7,154.58,126.19}; // critical temperatures
    double pc[4]={7.386e6,45.1e5,5.043e6,3.396e6}; // critical pressures
    double M[4]={44e-3,70e-3,32e-3,28e-3}; // molar masses
    double A[4][4];
    double gam[4];
    int i,j;
    for (i=0;i<x_size;i++) {
        y = y + x[i];
        k[i]=gas_thermal_conductivity[i];
    }
    for (i=0;i<x_size;i++) {
        xi[i] = x[i]/y; // normalize molar fractions
        gam[i]=210*pow(Tc[i]*pow(M[i],3)/pow(pc[i],4),1.0/6.0);
    }
    for (i=0;i<x_size;i++) {
        for (j=0;j<x_size;j++) {
            y=gam[j]*(exp(0.0464*T/Tc[i])-exp(-0.2412*T/Tc[i]))/\
                gam[i]/(exp(0.0464*T/Tc[j])-exp(-0.2412*T/Tc[j]));
            A[i][j]=pow(1+pow(y,0.5)*pow(M[i]/M[j],0.25),2)/\
                pow(8*(1+M[i]/M[j]),0.5);
        }
    }
    for (i=0;i<x_size;i++) {
        y=0;
        for (j=0;j<x_size;j++) {
            y=y+xi[j]*A[i][j];
        }
        kgasmix = kgasmix + xi[i]*k[i]/y;
    }
    outputs[0] = kgasmix;
}
'''
## Lindsay-Bromley model
#
#  see [link](http://dx.doi.org/10.1021/ie50488a017)
lindsayBromleyCode=r'''
#include "modena.h"
#include "math.h"
#include "stdio.h"

void gasMixtureConductivity
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    double kgasmix=0; // gas mixture conductivity
    double Rg=8.314; // gas constant
    double y=0;
    double xi[4]; // normalized molar fractions
    double k[4]; // thermal conductivities
    double Tb[4]={194.75,322.4,90.19,77.36}; // boiling point temperatures
    double cp[4]; // thermal capacities at constant pressure
    double M[4]={44e-3,70e-3,32e-3,28e-3}; // molar masses
    double S[4]; // Sutherland constants
    double gam[4]; // heat capacity ratio
    double cv[4]; // thermal capacities at constant volume
    double A[4][4];
    double t=T/1e3;
    int i,j;
    cp[0]=24.99735+55.18696*t-33.69137*t*t+7.948387*t*t*t-0.136638/t/t; // CO2
    cp[1]=-25.6132057+226.4176882*t+574.2688767*t*t-670.5517907*t*t*t+\
        0.6765321/t/t; // CyP
    cp[2]=31.32234-20.23531*t+57.86644*t*t-36.50624*t*t*t-0.007374/t/t; // O2
    cp[3]=28.98641+1.853978*t-9.647459*t*t+16.63537*t*t*t+0.000117/t/t; // N2
    for (i=0;i<x_size;i++) {
        y = y + x[i];
        k[i]=gas_thermal_conductivity[i];
    }
    for (i=0;i<x_size;i++) {
        // normalize molar fractions
        xi[i] = x[i]/y;
        S[i] = 1.5*Tb[i];
        cv[i] = cp[i] - Rg;
        gam[i] = cp[i]/cv[i];
    }
    for (i=0;i<x_size;i++) {
        for (j=0;j<x_size;j++) {
            y=k[i]/k[j]*cp[j]/cp[i]*(9-5/gam[j])/(9-5/gam[i]);
            A[i][j]=0.25*pow(1+pow(y*pow(M[j]/M[i],0.75)*(T+S[i])/\
                (T+S[j]),0.5),2)*(T+sqrt(S[i]*S[j]))/(T+S[j]);
        }
    }
    for (i=0;i<x_size;i++) {
        y=0;
        for (j=0;j<x_size;j++) {
            y=y+xi[j]*A[i][j];
        }
        kgasmix = kgasmix + xi[i]*k[i]/y;
    }
    outputs[0] = kgasmix;
}
'''

## Surrogate function for thermal conductivity of blowing agents.
#
# Thermal conductivity of blowing agents is a function of temperature.
f_gasMixtureConductivity = CFunction(
    Ccode=dohrnCode,
    # These are global bounds for the function
    inputs={
        'T': {'min': 273, 'max': 550},
        'x': {'index': gasConductivity.species, 'min': 0, 'max': 1},
        'gas_thermal_conductivity': {'index': gasConductivity.species,
            'min': 0, 'max': 1},
    },
    outputs={
        'gasMixtureConductivity': {'min': 0, 'max': +9e99, 'argPos': 0},
    },
    parameters={
        'param0': {'min': -9e99, 'max': +9e99, 'argPos': 0},
        'param1': {'min': -9e99, 'max': +9e99, 'argPos': 1},
        'param2': {'min': -9e99, 'max': +9e99, 'argPos': 2},
    },
)

## Surrogate model for thermal conductivity of mixture of blowing agents
#
# Forward mapping model is used.
m_gasMixtureConductivity = ForwardMappingModel(
    _id='gasMixtureConductivity',
    surrogateFunction=f_gasMixtureConductivity,
    substituteModels=[gasConductivity.m_CO2_thermal_conductivity,\
                      gasConductivity.m_CyP_thermal_conductivity,\
                      gasConductivity.m_O2_thermal_conductivity,\
                      gasConductivity.m_N2_thermal_conductivity],
    parameters=[1, 1, 1, 1],
)
