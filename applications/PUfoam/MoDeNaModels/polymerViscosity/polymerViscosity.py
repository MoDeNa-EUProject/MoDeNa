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
@file      Surrogate function and model definitions for polymer viscosity model.
@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

import os
import json
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

## Surrogate function for polymer viscosity.
#
# Polymer viscosity is a function of temperature and conversion.
f_polymerViscosity = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"

void viscosity_SM
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double Aeta = parameters[0];
    const double Eeta = parameters[1];
    const double AA = parameters[2];
    const double B = parameters[3];
    const double Xg = parameters[4];

    const double Rg = 8.31446218;

    if (X<Xg) {
        outputs[0] = Aeta*exp(Eeta/(Rg*T))*pow(Xg/(Xg-X),AA+B*X);
    } else {
        outputs[0] = 1e10;
    }
}
''',
    # These are global bounds for the function
    inputs={
        'T': {'min': 200, 'max': 550 },
        'X': {'min': 0, 'max': 1 },
    },
    outputs={
        'mu': {'min': 0, 'max': +9e99, 'argPos': 0},
    },
    parameters={
        'param1': {'min': -1e9, 'max': 1e9, 'argPos': 0},
        'param2': {'min': -1e9, 'max': 1e9, 'argPos': 1},
        'param3': {'min': -1e9, 'max': 1e9, 'argPos': 2},
        'param4': {'min': -1e9, 'max': 1e9, 'argPos': 3},
        'param5': {'min': -1e9, 'max': 1e9, 'argPos': 4},
    },
)
with open(os.getcwd()+'/inputs/unifiedInput.json') as jsonfile:
    inputs=json.load(jsonfile)
    X_gel=inputs['kinetics']['gelPoint']

## [literature data](http://dx.doi.org/10.1002/aic.690280213)
par = [4.1e-8, 38.3e3, 4.0, -2.0, X_gel]

## [literature data](http://dx.doi.org/10.1002/aic.690280213)
par2 = [10.3e-8, 41.3e3, 1.5, 1.0, X_gel]

## [literature data][1]
## [1]: http://dx.doi.org/10.1002/(SICI)1097-4628(19961017)62:3<567::AID-APP14>3.0.CO;2-W
par3 = [3.32e-8, 42.9e3, 2.32, 1.4, X_gel]

## based on [literature data](http://dx.doi.org/10.1002/aic.690280213), but
## gel point changed to 0.5 (Baser and Khakhar)
par4 = [4.1e-8, 38.3e3, 4.0, -2.0, X_gel]

## [literature data](http://doi.wiley.com/10.1002/aic.12002)
par5 = [3.1e0, 2.24e3, 3.5, -2.0, X_gel]

## [literature data](http://doi.wiley.com/10.1002/pen.760311605)
par6 = [1.6e-7, 44.9e3, 1.29, 1.86, X_gel]

## Surrogate model for polymer viscosity
#
# Forward mapping model is used.
m_polymerViscosity = ForwardMappingModel(
    _id='polymerViscosity',
    surrogateFunction=f_polymerViscosity,
    substituteModels=[],
    parameters=par6,
)
