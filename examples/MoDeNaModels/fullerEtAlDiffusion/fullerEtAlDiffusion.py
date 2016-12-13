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
@file      Implementation of diffuson coefficient model.
@author    Henrik Rusche
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   twoTank
"""

from modena import CFunction, IndexSet, ForwardMappingModel
import modena.Strategy as Strategy

species = IndexSet(
    name= 'species',
    names= [ 'H2O', 'N2', 'SO2' ]
)


f = CFunction(
    inputs={
        'T': { 'min': 0, 'max': 9e99 },
        'p': { 'min': 0, 'max': 9e99 },
    },
    outputs={
        'D[A]': { 'min': 0, 'max': 9e99, 'argPos': 0 },
    },
    parameters={
        'W[A]': { 'min': 0, 'max': 9e99, 'argPos': 0 },
        'V[A]': { 'min': 0, 'max': 9e99, 'argPos': 1 },
        'W[B]': { 'min': 0, 'max': 9e99, 'argPos': 2 },
        'V[B]': { 'min': 0, 'max': 9e99, 'argPos': 3 },
    },
    indices={
        'A': species,
        'B': species,
    },
    Ccode= '''
#include "modena.h"
#include "math.h"

void fullerEtAlDiffusion
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double WA = parameters[0];
    const double VA = parameters[1];
    const double WB = parameters[2];
    const double VB = parameters[3];

    outputs[0] = 1.011e-4*pow(T, 1.75)*pow(1.0/WA + 1.0/WB, 1.0/2.0);
    outputs[0] /= p*(pow(pow(VA, 1.0/3.0) + pow(VB, 1.0/3.0), 2.0));
}
''',
)

m = ForwardMappingModel(
    _id= 'fullerEtAlDiffusion[A=H2O,B=N2]',
    surrogateFunction= f,
    substituteModels= [ ],
    parameters= [ 16, 9.44, 14, 11.38 ],
)
