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
Surrogate function and model definitions for the density of the reaction mixturemode
@author    Erik Laurini
@author    Mohsen Karimi
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
import modena
from modena import CFunction, IndexSet, \
    ForwardMappingModel, BackwardMappingModel, SurrogateModel
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
from jinja2 import Template

## Create terminal for colour output
term = Terminal()
## Surrogate function for density of the reaction mixture.
#
# Density is a function of temperature and conversion of gelling reaction.
f = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"
void densityreactionmixture
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double a_L0 = parameters[0];
    const double b_L0 = parameters[1];
    const double a_P0 = parameters[2];
    const double b_P0 = parameters[3];

    outputs[0] = ((a_L0*T + b_L0) + ((a_P0*T + b_P0) - (a_L0*T + b_L0))*XOH);
}
''',
    # These are global bounds for the function
    inputs={
        'T': { 'min': 273, 'max': 550},
        'XOH': { 'min': 0, 'max': 1},
    },
    outputs={
        'density_polymer': { 'min': 0, 'max': 8000, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -9e99, 'max': +9e99, 'argPos': 0 },
        'param1': { 'min': 0.0, 'max': +9e99, 'argPos': 1 },
        'param2': { 'min': -9e99, 'max': +9e99, 'argPos': 2 },
        'param3': { 'min': 0.0, 'max': +9e99, 'argPos': 3 },
    },
)
## Surrogate model for density_reaction_mixture
#
# Forward mapping model is used.
m = ForwardMappingModel(
    _id='density_reaction_mixture',
    surrogateFunction=f,
    substituteModels= [ ],
    parameters=[-0.0006, 1287.8, -0.000048, 992.8],
)
