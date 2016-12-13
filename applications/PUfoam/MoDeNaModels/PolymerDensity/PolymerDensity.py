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
@todo      Document
@author    Jonas Mairhofer
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

import os
import modena
from modena import ForwardMappingModel,BackwardMappingModel,SurrogateModel,CFunction,IndexSet,ModenaFireTask
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
from jinja2 import Template

# Create terminal for colour output
term = Terminal()


blowing_agents = IndexSet(
    name= 'blowing_agents',
    names= [ 'air', 'CO2']
)

monomers = IndexSet(
    name = 'monomers',
    names = ['PU', 'THF', 'hexane', 'surfactant']
)


# ********************************* Class ********************************** #
@explicit_serialize
class DensityExactSim(ModenaFireTask):
    """
    A FireTask that starts a microscopic code and updates the database.
    """

    def task(self, fw_spec):
        # Generate input fileblock
        self.generate_inputfile()

        # Execute detailed model
        ret = os.system(os.path.dirname(os.path.abspath(__file__))+'/src/PCSAFT_Density')

        # Check framework for errors
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

        #create output file for detailed code
        with open('out.txt', 'w+') as FILE:
            pass


    def analyse_output(self):
        """Method analysing the output of the file.
        @TODO consider adding check for empty file
        """
        with open('out.txt', 'r') as FILE:
            self['point']['rho'] = float(FILE.readline())


f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void surroDensity
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

    const double expo = 1.0 + (1.0 - T/P2);
    const double pwr  = pow(P1,expo);

    outputs[0] = P0 / pwr;
    //outputs[0] = P0 + T*P1 + P2*T*T;
}
''',
    # These are global bounds for the function
    inputs={
        'T': { 'min': 270.0, 'max': 550.0},
    },
    outputs={
        'rho': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -1E10, 'max': 1E10+2, 'argPos': 0 },    #check if boundaries are reasonable!!!
        'param1': { 'min': -1E10, 'max': 1E10+2, 'argPos': 1 },
        'param2': { 'min': -1E10, 'max': 1E10+2, 'argPos': 2 },
    },
    species={
        'A' : blowing_agents,
        'B' : monomers
    }
)

m = BackwardMappingModel(
    _id= 'PolymerDensity[A=AIR,B=PU]',
    surrogateFunction= f,
    exactTask= DensityExactSim(),
    substituteModels= [ ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [270.0, 300.0, 350.0, 400, 450, 500, 550],
                 },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 10.0,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)
