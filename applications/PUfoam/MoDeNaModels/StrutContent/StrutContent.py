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
Surrogate function and model definitions for the strut content
@author    Pavel Ferkl
@author    Mohsen Karimi
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
import modena
from modena import CFunction, ForwardMappingModel
## Surrogate function for strut content.
#
# Strut content is a assumed to be a function of foam density.
f = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"
#define MAX(a,b) ((a) > (b) ? a : b)
void strutContent
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double a = parameters[0];
    const double b = parameters[1];
    const double c = parameters[2];
    const double pi = 3.14159;

    outputs[0] = MAX(atan(a*rho + b)/pi*2,0.0);
}
''',
    # These are global bounds for the function
    inputs={
        'rho': { 'min': 0, 'max': 1e5},
    },
    outputs={
        'fs': { 'min': 0, 'max': 1, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -9e99, 'max': +9e99, 'argPos': 0 },
        'param1': { 'min': -9e99, 'max': +9e99, 'argPos': 1 },
        'param2': { 'min': -9e99, 'max': +9e99, 'argPos': 2 },
    },
)
## Surrogate model for density_reaction_mixture
#
# Forward mapping model is used.
m = ForwardMappingModel(
    _id='strutContent',
    surrogateFunction=f,
    substituteModels= [ ],
    parameters=[0.06115509, -0.72513392,  1.],
)
