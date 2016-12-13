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
Surrogate function and model definitions for thermal conductivity of
polyurethane

@author    Erik Laurini
@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
import modena
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, CFunction
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
from jinja2 import Template

## Create terminal for colour output
term = Terminal()

## Surrogate function for thermal conductivity of polyurethane.
#
# Thermal conductivity of polyurethane is a function of temperature.
f_polymer_thermal_conductivity = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"

void thermal_conductivity
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double a = parameters[0];
    const double b = parameters[1];

    outputs[0] = (a*T)+b;
}
''',
    # These are global bounds for the function
    inputs={
        'T': {'min': 273, 'max': 550},
    },
    outputs={
        'polymer_thermal_conductivity': {
            'min': 0, 'max': +9e99, 'argPos': 0
        },
    },
    parameters={
        'param0': {'min': -9e99, 'max': +9e99, 'argPos': 0},
        'param1': {'min': -9e99, 'max': +9e99, 'argPos': 1},
    },
)

## Surrogate model for thermal conductivity of polyurethane
#
# Forward mapping model is used.
m_polymer_thermal_conductivity = ForwardMappingModel(
    _id='polymer_thermal_conductivity',
    surrogateFunction=f_polymer_thermal_conductivity,
    substituteModels=[],
    parameters=[0.198e-3, 151.08e-3],
    inputs={
        'T': {'min': 273, 'max': 450},
    },
    outputs={
        'polymer_thermal_conductivity': {'min': 0, 'max': +9e99},
    },
)
