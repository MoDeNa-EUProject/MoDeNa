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
@file      Backward mapping FireTask for rheology model.
@author    Christos Mitrias
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

import os
import modena
import SurfaceTension, polymerViscosity
import Rheology
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, CFunction, ModenaFireTask
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from jinja2 import Template

f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"
#include <stdio.h>

void Arrhenius
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double A_mu = 0.0387; // could be loaded from some library for consistence
    const double E_mu = 10000;
    const double R_rh = 8.314;

    double f_t, mu_ap;

    f_t = A_mu * exp(E_mu / R_rh / T );

    mu_ap = mu_car * f_t;
    // printf("apparent viscosity %f", mu_ap);
    outputs[0] = mu_ap;
}
''',
   inputs={
       'T': {'min': 0, 'max': 9e99 },
       'shear': {'min': 0, 'max': 9e99 },
       'X': {'min': 0, 'max': 1 },
       'm0' : {'min': 0, 'max' : 9e99},
       'm1' : {'min': 0, 'max' : 9e99},
       'mu' : {'min': 0, 'max' : 9e99},
       'ST' : {'min': 0, 'max' : 9e99},
       'mu_car' : {'min': 0, 'max' : 9e99},
   },
   outputs={
       'mu_ap': { 'min': 0, 'max': 9e99, 'argPos': 0 },
   },
   parameters={
       'Rgas' : {'min': 8.31, 'max': 8.32, 'argPos': 0},
   }
)


m = ForwardMappingModel(
    _id='Rheology_Arrhenius',
    surrogateFunction=f,
    substituteModels=[ Rheology.m ],
    parameters=[ 8.314 ],
)
