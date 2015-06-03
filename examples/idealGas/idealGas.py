'''

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014 MoDeNa Consortium, All rights reserved.

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

Description
    Python library of FireTasks

Authors
    Henrik Rusche

Contributors
'''

import os
import modena
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, CFunction
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal

# Create terminal for colour output
term = Terminal()


__author__ = 'Henrik Rusche'
__copyright__ = 'Copyright 2014, MoDeNa Project'
__version__ = '0.2'
__maintainer__ = 'Henrik Rusche'
__email__ = 'h.rusche@wikki.co.uk.'
__date__ = 'Sep 4, 2014'


f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void idealGas
(
    const double* parameters,
    const double* inherited_inputs,
    const double* inputs,
    double *outputs
)
{
    const double p0 = inputs[0];
    const double T0 = inputs[1];

    const double R = parameters[0];

    outputs[0] = p0/R/T0;
}
''',
    # These are global bounds for the function
    inputs={
        'p0': { 'min': 0, 'max': 9e99, 'argPos': 0 },
        'T0': { 'min': 0, 'max': 9e99, 'argPos': 1 },
    },
    outputs={
        'rho0': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'R': { 'min': 0.0, 'max': 9e99, 'argPos': 0 }
    },
)

m1 = ForwardMappingModel(
    _id= 'idealGas',
    #exactTask= Strategy.InitialDataPoints(),
    surrogateFunction= f,
    substituteModels= [ ],
    parameters= [ 287.0 ],
    inputs={
        'p0': { 'min': 0, 'max': 9e99 },
        'T0': { 'min': 0, 'max': 9e99 },
    },
    outputs={
        'rho0': {'min': 0, 'max': 9e99 },
    },
    initialisationStrategy= Strategy.InitialData(
        initialData={
        'p0' : [1,2,3,4,5],
        'T0' : [296,297,298,299,300],
        'rho0' : [0.000011771353234, 0.00002346343809758, 0.00003507705259219, 0.00004661298404, 0.0000580720092915],
        },
    ),
)

m2 = BackwardMappingModel(
    _id= 'idealGasBackward',
    exactTask= Strategy.InitialDataPoints(),
    surrogateFunction= f,
    substituteModels= [ ],
    parameters= [ 287.0 ],
    inputs={
        'p0': { 'min': 0, 'max': 9e99 },
        'T0': { 'min': 0, 'max': 9e99 },
    },
    outputs={
        'rho0': {'min': 0, 'max': 9e99 },
    },
    initialisationStrategy= Strategy.InitialData(
        initialData={
        'p0' : [1,2,3,4,5],
        'T0' : [296,297,298,299,300],
        'rho0' : [0.000011771353234, 0.00002346343809758, 0.00003507705259219, 0.00004661298404, 0.0000580720092915],
        },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.1,
        maxError= 0.05,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)


