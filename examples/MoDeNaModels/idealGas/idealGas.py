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
@file      Implementation of ideal gas model.
@author    Henrik Rusche
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   twoTank
"""

import os
import modena
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, CFunction, ModenaFireTask
import modena.Strategy as Strategy
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize


f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void idealGas
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double R = parameters[0];

    outputs[0] = p0/R/T0;
}
''',
    # These are global bounds for the function
    inputs={
        'p0': { 'min': 0, 'max': 9e99 },
        'T0': { 'min': 0, 'max': 9e99 },
    },
    outputs={
        'rho0': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'R': { 'min': 0.0, 'max': 9e99, 'argPos': 0 }
    },
)

m = ForwardMappingModel(
    _id= 'idealGas',
    surrogateFunction= f,
    substituteModels= [ ],
    parameters= [ 287.0 ],
)
