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
@namespace  FoamElasticModulus.ElasticModulus
@ingroup    mod_elasticmodulus
@brief      Surrogate Function, Surrogate Model templates and Model Recipe
@author
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@details

Detailed


description


of


the


module
"""

import os
from os import getcwd, system
from os.path import abspath, dirname, join
from jinja2 import Template
import subprocess

from modena import ForwardMappingModel,BackwardMappingModel, CFunction, ModenaFireTask, Strategy


from fireworks.utilities.fw_utilities import explicit_serialize

# ----------------------- Convenience variables ----------------------------- #
MODULE_DIR = dirname(abspath(__file__))
EXECUTION_DIR = getcwd()

TEMPLATE_FoamElasticModulus = Template(\
"""
{{ s['point']['Mu'] }}
{{ s['point']['Sigma'] }}
{{ s['point']['strut_content'] }}
{{ s['point']['rho'] }}
""", trim_blocks=True, lstrip_blocks=True)

# ----------------------- Exact simulation wrapper -------------------------- #
@explicit_serialize
class MechanicalPropertiesExactSim(ModenaFireTask):
    """
    """

    INPUT_FILE = "Inputs.in"
    OUTPUT_FILE = "Outputs.in"

    def task(self, fw_spec):

        # Write input for detailed model
        TEMPLATE_FoamElasticModulus.stream(s=self).dump(self.INPUT_FILE)

        # Execute detailed model
        run_command = join(MODULE_DIR, "AbaqusSimulation.py")

        ret = subprocess.check_call(run_command, shell=True)
        # This call enables backward mapping capabilities (not needed in this example)
        self.handleReturnCode(ret)

        # Analyse output
        with open(self.OUTPUT_FILE, 'r') as FILE:
            self['point']['E'] = float(FILE.readline())



# ----------------------- Surrogate Functions ------------------------------- #
f_backward = CFunction(
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


    const double C1 = parameters[0];
    const double C2 = parameters[1];

    const double E_PU = 2400;
    const double rho_PU = 1200;

    outputs[0] = E_PU*(C1*pow(strut_content*rho/rho_PU,2) + C2*(1-strut_content)*(rho/rho_PU));
}
''',
    # These are global bounds for the function
    inputs={
        'rho': { 'min': 0, 'max': 1e5},
        'strut_content': { 'min': 0, 'max': 1e5},
        'Mu': { 'min': 0, 'max': 1e5},
        'Sigma': { 'min': 0, 'max': 1e5},
    },
    outputs={
        'E': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'C1': { 'min': 0.0, 'max': 10.0, 'argPos': 0 },
        'C2': { 'min': 0.0, 'max': 10.0, 'argPos': 1 },
    },
)


# Strut content is a assumed to be a function of foam density.
f_ludwigshafen = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"
void ElasticModulus
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double C1 = parameters[0];
    const double C2 = parameters[1];
    const double E_PU = 2400;
    const double rho_PU = 1200;

    outputs[0] = E_PU*(C1*pow(strut_content*rho/rho_PU,2) + C2*(1-strut_content)*(rho/rho_PU));
}
''',
    # These are global bounds for the function
    inputs={
        'rho': { 'min': 0, 'max': 1e5},
        'strut_content': { 'min': 0, 'max': 1e5},
        'Mu': { 'min': 0, 'max': 1e5},
        'Sigma': { 'min': 0, 'max': 1e5},
    },
    outputs={
        'E': { 'min': 0, 'max': 1, 'argPos': 0 },
    },
    parameters={
        'C1': { 'min': -9e99, 'max': +9e99, 'argPos': 0 },
        'C2': { 'min': -9e99, 'max': +9e99, 'argPos': 1 },
    },
)


# -------------------------- Surrogate Models ------------------------------- #
m_forward = ForwardMappingModel(
    _id='FoamElasticModulus',
    surrogateFunction=f_ludwigshafen,
    substituteModels= [ ],
    parameters=[0, 1],# TODO: Use real values from Ludwigshafen presentation
)


m_backward = BackwardMappingModel(
    _id= 'FoamElasticModulus_TEST',
    surrogateFunction= f_backward,
    exactTask=MechanicalPropertiesExactSim(),
    substituteModels= [],
    initialisationStrategy= Strategy.InitialData(
       initialData={
        'rho': [30.2,60,120,240,30.2,60,120,240,30.2,30.2,30.2,30.2,30.2,30.2,30.2],
        'strut_content': [0.4,0.4,0.4,0.4,0.8,0.8,0.8,0.8,0.85,0.85,0.85,0.85,0.85,0.85,0.85],
        'Mu': [0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39,0.39],
        'Sigma': [0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.01,0.05,0.1,0.15,0.2,0.25,0.30],
        'E': [8.16,22.85,52.89,124.73,6.53,13.06,30.36,65.63,4.97,5.48,5.42,4.52,4.61,4.42,4.79],
        }
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 9e99,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)

