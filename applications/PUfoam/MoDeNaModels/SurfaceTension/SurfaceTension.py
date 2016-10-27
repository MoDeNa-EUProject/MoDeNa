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
This is the Surface Tension python module. Basically, it contains the following:

The FireTask which controls the call of the detailed model. This detailed model is called
at the very beginning of the simulation in order to generate initial data points
which can be used to fit the parameters of the surrogate model and during a running simulation
as soon as the Surface Tension model is called with input parameters which lie outside the range
the parameters of the surrogate model was so far fitted for. This FireTask is stored in the class
"SurfaceTensionExactSim" and a more detailed description of the detailed model can be found
in the description of this class.

Furthermore, this module contains the code of the surrogate model function as well as the
definitions of its input and output values and its fittable parameters. Care should be
taken to set reasonable bounds for these variables.

Also, this module contains the backward mapping model. This model consits of the
surrogate model function, an initialisation strategy, the out of bounds strategy and the
parameter fitting strategy. The initialisation strategy defines the initial data points where the
detailed model will be evaluated at simulation start for an initial fit of the surrogate model parameters.
The out of bounds strategy determines, how many new points and where to place these new
points, once the Surface Tension model is called for input values outside of the
fitted range. The parameter fitting strategy defines tolerances and maximal iterations
which are passed to the numerical solver which performs the actual fitting of the
surrogate model parameters.


@author    Jonas Mairhofer
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

import os
from modena import *
from fireworks.utilities.fw_utilities import explicit_serialize
from jinja2 import Template


blowing_agents = IndexSet(
    name= 'blowing_agents',
    names= [ 'AIR', 'CO2']
)

monomers = IndexSet(
    name = 'monomers',
    names = ['PU', 'THF', 'HEXANE', 'surfactant']
)



# ********************************* Class ********************************** #
@explicit_serialize
class SurfaceTensionExactSim(ModenaFireTask):
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

    def task(self, fw_spec):

        # Write input for detailed model
        self.generate_inputfile()

        # Execute detailed model
        run_command = os.path.dirname(os.path.abspath(__file__))+'/src/PCSAFT_SurfaceTension -snes_monitor_short -ksp_monitor_short \
                   -nx 800 -rc 9.0 -box 300 -erel 1e-08 -init_pert 0 \
                   -snes_type newtonls  -snes_converged_reason  \
                   -snes_atol 1e-07 -snes_rtol 1e-07 -snes_stol 1e-07 -snes_max_it 15 \
                   -ksp_max_it 10 -ksp_gmres_restart 50 \
                   -snes_linesearch_type l2 -snes_linesearch_damping 0.3 -snes_linesearch_monitor \
                   -snes_max_fail 1 -snes_max_linear_solve_fail 100 \
                   -ksp_gmres_cgs_refinement_type refine_always \
                   -jac 0 -pc_type none > log'


        ret = os.system(run_command)
        # This call enables backward mapping capabilities (not needed in this example)
        self.handleReturnCode(ret)

        # Analyse output
        self.analyse_output()

    def generate_inputfile(self):
        """Method generating a input file using the Jinja2 template engine."""
        Template("""
            {#
                 Write inputs to the template, one per line.
            #}
            {% for k,v in s['point'].iteritems() %}
                 {{ v }}
            {% endfor %}
            {#
                 The number of species, one integer.
            #}
                  {{ s['indices'].__len__() }}
            {#
                 Write the species (lower case) one per line.
            #}
            {% for k,v in s['indices'].iteritems() %}
                   {{ v.lower() }}
            {% endfor %}
            {#
                    Set initial feed molar fractions to zero.
            #}
            {% for k,v in s['indices'].iteritems() %}
                {{ 0.0 }}
            {% endfor %}
            """, trim_blocks=True,
               lstrip_blocks=True).stream(s=self).dump('in.txt')

        with open('out.txt','w+') as FILE:
            pass

    def analyse_output(self):
        """Method        'A' : blowing_agents,
        'B' : monomers analysing the output of the file.
             @TODO consider adding check for empty file
        """
        with open('out.txt', 'r') as FILE:
            self['point']['ST'] = float(FILE.readline())




f = CFunction(
    Ccode= r'''
#include "modena.h"
#include "math.h"

void surroSurfaceTension
(
const modena_model_t* model,
const double* inputs,
double *outputs
)
{
{% block variables %}{% endblock %}

const double P0 = parameters[0];
const double P1 = parameters[1];
const double P2 = parameters[2];

outputs[0] = P0 + T*P1 + P2*T*T;
}
''',
    # These are global bounds for the function
    inputs={
        'T': { 'min': 270.0, 'max': 550.0 },        #check if boundaries reasonable, from this range, the random values for the DOE are chosen!
    },
    outputs={
        'ST': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -1E10, 'max': 1E10, 'argPos': 0 },    #check if boundaries are reasonable!!!
        'param1': { 'min': -1E10, 'max': 1E10, 'argPos': 1 },
        'param2': { 'min': -1E10, 'max': 1E10, 'argPos': 2 },
    },
    species = {
        'A' : blowing_agents,
        'B' : monomers,
    }
)

m = BackwardMappingModel(
    _id= 'SurfaceTension[A=AIR,B=THF]',
    surrogateFunction= f,
    exactTask= SurfaceTensionExactSim(),
    substituteModels= [ ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [270.0, 290.0, 330.0],
        },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 1e-2,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)

m2 = BackwardMappingModel(
    _id= 'SurfaceTension[A=AIR,B=PU]',
    surrogateFunction= f,
    exactTask= SurfaceTensionExactSim(),
    substituteModels= [ ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [340.0, 350.0, 360.0, 371.0, 380.0],
        },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 1e-0,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)

