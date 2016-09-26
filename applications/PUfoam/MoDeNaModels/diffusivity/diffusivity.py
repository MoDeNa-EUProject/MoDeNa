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
Surrogate function and model definitions for diffusivity of blowing agents in
polymer.

@author    Erik Laurini
@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
from modena import *
from fireworks.utilities.fw_utilities import explicit_serialize
from jinja2 import Template


## List of components, for which surrogate model is provided
species = IndexSet(
    name= 'diffusivity_pol_species',
    names= [ 'CO2', 'CyP', 'N2', 'O2' ]
)
## Surrogate function for diffusivity of blowing agents in polymer.
#
# Diffusivity is a function of temperature.
f_diffusivity = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"

void diffusivityPol
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double a = parameters[0];
    const double b = parameters[1];

    outputs[0] = a*exp(-(b*(1/T)));
}
''',
    # These are global bounds for the function
    inputs={
        'T': {'min': 273, 'max': 550},
    },
    outputs={
        'diffusivity': {'min': 0, 'max': +9e99, 'argPos': 0},
    },
    parameters={
        'param0[A]': {'min': 0.0, 'max': +9e99, 'argPos': 0},
        'param1[A]': {'min': 0.0, 'max': +9e99, 'argPos': 1},
    },
    indices={
        'A': species,
    },
)
## Surrogate model for diffusivity
#
# Forward mapping model is used.
m_CO2_diffusivity = ForwardMappingModel(
    _id='diffusivityPol[A=CO2]',
    surrogateFunction=f_diffusivity,
    substituteModels=[],
    parameters=[0.00123, 6156],
)
## Surrogate model for diffusivity
#
# Forward mapping model is used.
m_CyP_diffusivity = ForwardMappingModel(
    _id='diffusivityPol[A=CyP]',
    surrogateFunction=f_diffusivity,
    substituteModels=[],
    parameters=[1.7e-7, 4236],
)
## Surrogate model for diffusivity
#
# Forward mapping model is used.
m_N2_diffusivity = ForwardMappingModel(
    _id='diffusivityPol[A=N2]',
    surrogateFunction=f_diffusivity,
    substituteModels=[],
    parameters=[0.003235, 6927],
)
## Surrogate model for diffusivity
#
# Forward mapping model is used.
m_O2_diffusivity = ForwardMappingModel(
    _id='diffusivityPol[A=O2]',
    surrogateFunction=f_diffusivity,
    substituteModels=[],
    parameters=[0.00085, 6411],
)
