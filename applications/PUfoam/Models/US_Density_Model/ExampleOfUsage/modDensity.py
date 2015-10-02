
'''

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014-2015 MoDeNa Consortium, All rights reserved.

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
from modena import ForwardMappingModel,BackwardMappingModel,SurrogateModel,CFunction,IndexSet
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
from jinja2 import Template

# Create terminal for colour output
term = Terminal()


__author__ = 'Henrik Rusche'
__copyright__ = 'Copyright 2014, MoDeNa Project'
__version__ = '0.2'
__maintainer__ = 'Henrik Rusche'
__email__ = 'h.rusche@wikki.co.uk.'
__date__ = 'Sep 4, 2014'


# ********************************* Class ********************************** #
@explicit_serialize
class DensityExactSim(FireTaskBase):
    """
    A FireTask that starts a microscopic code and updates the database.
    """

    def run_task(self, fw_spec):
        print(
            term.yellow
          + "Performing exact simulation (microscopic code recipe)"
          + term.normal
        )

        # Write input for detailed model
        ff = open('in.txt', 'w')
        Tstr = str(self['point']['T'])
        ff.write('%s \n' %(Tstr))


        ##TODO INPUT SHOULD COME FROM IndexSet

        ff.write('2 \n')       #number of components in system
        ff.write('air \n')     #component 1
        ff.write('pu \n')      #component 2
        ff.write('0. \n')     #molar feed (initial) composition component 1
        ff.write('0. \n')     #molar feed (initial) composition component 2
        ff.close()

        #create output file for detailed code
        fff = open('out.txt', 'w+')
        fff.close()

        # Execute detailed model
        os.system('../src/PCSAFT_Density')

        # Analyse output
        f = open('out.txt', 'r')
        self['point']['rho'] = float(f.readline())
        f.close()

        return FWAction(mod_spec=[{'_push': self['point']}])






f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void surroDensity
(
    const double* parameters,
    const double* inherited_inputs,
    const double* inputs,
    double *outputs
)
{

    const double T = inputs[0];

    const double P0 = parameters[0];
    const double P1 = parameters[1];
    const double P2 = parameters[2];

   // const double expo = 1.0 + (1.0 - T/P2);
   // const double pwr  = pow(P1,expo);

   // outputs[0] = P0 / pwr;


    outputs[0] = P0 + T*P1 + P2*T*T;
}
''',
    # These are global bounds for the function
    inputs={
        'T': { 'min': 270.0, 'max': 300.0, 'argPos': 0 },        #check if boundaries reasonable, from this range, the random values for the DOE are chosen!
    },
    outputs={
        'rho': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -1E10, 'max': 1E10, 'argPos': 0 },    #check if boundaries are reasonable!!!
        'param1': { 'min': -1E10, 'max': 1E10, 'argPos': 1 },
        'param2': { 'min': -1E10, 'max': 1E10, 'argPos': 2 },
    },
)

m = BackwardMappingModel(
    _id= 'Density',
    surrogateFunction= f,
    exactTask= DensityExactSim(),
    substituteModels= [ ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [270.0, 290.0, 300.0],
                 },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 300.0,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)

