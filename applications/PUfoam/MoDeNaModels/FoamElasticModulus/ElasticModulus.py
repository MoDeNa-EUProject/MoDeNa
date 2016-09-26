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

import os
from modena import ForwardMappingModel,BackwardMappingModel,SurrogateModel,CFunction,IndexSet#,ModenaFireTask
from fireworks.utilities.fw_utilities import explicit_serialize

'''
# ********************************* Class ********************************** #
@explicit_serialize
class MechanicalPropertiesExactSim(ModenaFireTask):
    """
    This FireTask controls the execution of the detailed model of the Surface Tension model.
    The detailed model is a density functional theory implementation based on PC-SAFT. A
    detailed description of this model can be found in Deliverable 1.3 on the MoDeNa website.
    In order to start the detailed model, the input values for the model are first written to the
    file "in.txt". The detailed model code picks them up from this file and performs the according
    calculation. Once it is done, the output value is written to the file "out.txt". This FireTask
    then reads in the calculated surface tension from "out.txt" and inserts this value into the
    database.
    """

    MODULE_DIR = os.path.dirname(os.path.abspath(__file__))#ectory containing this file
    INPUT_FILE = "Inputs.in"
    OUTPUT_FILE = "Outputs.in"

    def task(self, fw_spec):

        # Write input for detailed model
        self.generate_inputfile()

        # Execute detailed model
        # run_command = MODULE_DIR + "/python FoamElasticModulus.py"
        run_command = MODULE_DIR + "/python 'Hello World'"

        ret = os.system(run_command)
        # This call enables backward mapping capabilities (not needed in this example)
        self.handleReturnCode(ret)

        # Analyse output
        self.analyse_output()

    def generate_inputfile(self):
        """Method generating a input file using the Jinja2 template engine."""
        Template("""
{{ s['point']['Mu'] }}
{{ s['point']['Sigma'] }}
            """, trim_blocks=True,
               lstrip_blocks=True).stream(s=self).dump(self.INPUT_FILE)

    def analyse_output(self):
        """Method analysing the output of the file.
             @TODO consider adding check for empty file
        """
        with open(self.OUTPUT_FILE, 'r') as FILE:
            self['point']['E'] = float(FILE.readline())
'''
"""
# TODO: This function is not implemented properly
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
}
''',
    # These are global bounds for the function
    inputs={
        'Mu': { 'min': 0, 'max': 9e99 },
        'Var': { 'min': 0, 'max': 9e99 },
    },
    outputs={
        'E': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': 0.0, 'max': 10.0, 'argPos': 0 },
        'param1': { 'min': 0.0, 'max': 10.0, 'argPos': 1 },
    },
)
"""


# Strut content is a assumed to be a function of foam density.
f = CFunction(
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


# Forward mapping model is used.
m_forward = ForwardMappingModel(
    _id='FoamElasticModulus',
    surrogateFunction=f,
    substituteModels= [ ],
    parameters=[0, 1],# TODO: Use real values from Ludwigshafen presentation
)
"""
m_backward = BackwardMappingModel(
    _id= 'FoamElasticModulus',
    surrogateFunction= f,
    exactTask= ,
    substituteModels= [ ],
    initialisationStrategy= Strategy.EmptyInitialisationStrategy(),
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
"""
