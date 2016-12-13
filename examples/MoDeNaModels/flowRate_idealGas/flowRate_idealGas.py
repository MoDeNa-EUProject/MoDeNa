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
@file      Implementation of flow rate model with substituted ideal gas.
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
from blessings import Terminal
from jinja2 import Template
import idealGas


# ********************************* Class ********************************** #
@explicit_serialize
class FlowRateExactSim(ModenaFireTask):
    """
    A FireTask that starts a microscopic code and updates the database.
    """

    def task(self, fw_spec):
        # Write input

        # See http://jinja.pocoo.org/docs/dev/templates/
        Template('''
{{ s['point']['D'] }}
{{ s['point']['rho0'] }}
{{ s['point']['p0'] }}
{{ s['point']['p1Byp0'] }}
        '''.strip()).stream(s=self).dump('in.txt')

        # Execute the application
        # In this simple example, this call stands for a complex microscopic
        # code - such as full 3D CFD simulation.
        # Source code in src/flowRateExact.C
        ret = os.system(os.path.dirname(os.path.abspath(__file__))+'/src/flowRateExact')

        # This enables backward mapping capabilities (not needed in this example)
        self.handleReturnCode(ret)

        # Analyse output
        f = open('out.txt', 'r')
        self['point']['flowRate'] = float(f.readline())
        f.close()


f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void two_tank_flowRate
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double D = inputs[0];
    const double rho0 = inputs[1];
    const double p0 = inputs[2];
    const double p1 = p0*inputs[3];

    const double P0 = parameters[0];
    const double P1 = parameters[1];

    outputs[0] = M_PI*pow(D, 2.0)*P1*sqrt(P0*rho0*p0);
}
''',
    # These are global bounds for the function
    inputs={
        'D': { 'min': 0, 'max': 9e99 },
        'T0': { 'min': 0, 'max': 9e99 },
        'p0': { 'min': 0, 'max': 9e99 },
        'p1Byp0': { 'min': 0, 'max': 1.0 },
    },
    outputs={
        'flowRate': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': 0.0, 'max': 10.0, 'argPos': 0 },
        'param1': { 'min': 0.0, 'max': 10.0, 'argPos': 1 },
    },
)

m = BackwardMappingModel(
    _id= 'flowRate',
    surrogateFunction= f,
    exactTask= FlowRateExactSim(),
    substituteModels= [ idealGas.m ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'D': [0.01, 0.01, 0.01, 0.01],
            'T0': [300, 300, 300, 300],
            'p0': [2.8e5, 3.2e5, 2.8e5, 3.2e5],
            'p1Byp0': [0.03, 0.03, 0.04, 0.04],
        },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 0.05,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)
