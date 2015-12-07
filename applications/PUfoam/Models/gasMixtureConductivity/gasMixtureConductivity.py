'''@cond

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
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.
@endcond'''

"""
@file
Surrogate function and model definitions for thermal conductivity of mixture of
blowing agents.

@author    Pavel Ferkl
@copyright 2014-2015, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
import modena
from modena import CFunction, IndexSet, Workflow2, \
    ForwardMappingModel, BackwardMappingModel, SurrogateModel
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
import gasConductivity

## Create terminal for colour output
term = Terminal()

## Surrogate function for thermal conductivity of blowing agents.
#
# Thermal conductivity of blowing agents is a function of temperature.
f_gasMixtureConductivity = CFunction(
    Ccode=r'''
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
    int i;

    printf("temp = %g\\n", T);
    printf("x = ");
    for (i=0;i<x_size;i++) {
        printf("%g ", x[i]);
    }
    printf("\n");
    printf("k = ");
    for (i=0;i<x_size;i++) {
        printf("%g ", gas_thermal_conductivity[i]);
    }
    printf("\n");

    for (i=0;i<x_size;i++) {
        kgasmix = kgasmix + parameters[i]*gas_thermal_conductivity[i]*x[i];
    }
    outputs[0] = kgasmix;
}
''',
    # These are global bounds for the function
    inputs={
        'T': {'min': 273, 'max': 450},
        'x': {'index': gasConductivity.species, 'min': 0, 'max': 1},
        'gas_thermal_conductivity': {'index': gasConductivity.species, 'min': 0, 'max': 1},
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
                      gasConductivity.m_Air_thermal_conductivity,\
                      gasConductivity.m_CyP_thermal_conductivity],
    parameters=[1, 1, 1],
)
