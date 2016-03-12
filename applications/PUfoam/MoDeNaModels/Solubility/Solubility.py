'''@cond

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
@endcond'''

"""
@file
Python library of FireTasks
This is the Solubility python module. Basically, it contains the following:

The FireTask which controls the call of the detailed model. This detailed model is called
at the very beginning of the simulation in order to generate initial data points
which can be used to fit the parameters of the surrogate model and during a running simulation
as soon as the Solubility model is called with input parameters which lie outside the range
the parameters of the surrogate model was so far fitted for. This FireTask is stored in the class
"SolubilityExactSim" and a more detailed description of the detailed model can be found
in the description of this class.

Furthermore, this module contains the code of the surrogate model function as well as the
definitions of its input and output values and its fittable parameters. Care should be
taken to set reasonable bounds for these variables.

Also, this module contains the backward mapping model. This model consits of the
surrogate model function, an initialisation strategy, the out of bounds strategy and the
parameter fitting strategy. The initialisation strategy defines the initial data points where the
detailed model will be evaluated at simulation start for an initial fit of the surrogate model parameters.
The out of bounds strategy determines, how many new points and where to place these new
points, once the Solubility model is called for input values outside of the
fitted range. The parameter fitting strategy defines tolerances and maximal iterations
which are passed to the numerical solver which performs the actual fitting of the
surrogate model parameters.

@author    Jonas Mairhofer, Pavel Ferkl
@copyright 2014-2015, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

import os
import modena
from modena import ForwardMappingModel,BackwardMappingModel,SurrogateModel,\
    CFunction,IndexSet
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
from jinja2 import Template

# Create terminal for colour output
term = Terminal()

species = IndexSet(
    name= 'solubility_pol_species',
    names= [ 'Air', 'CO2', 'CyP' ]
)

@explicit_serialize
class SolubilityExactSim(FireTaskBase):
    """
    This FireTask controls the execution of the detailed model of the Solubility model.
    The detailed model uses the PC-SAFT equation of state. A
    detailed description of PC-SAFT model can be found in Deliverable 1.3 on the MoDeNa website.

    In order to start the detailed model, the input values for the model are first written to the
    file "in.txt". The detailed model code picks them up from this file and performs the according
    calculation. Once it is done, the output value is written to the file "out.txt". This FireTask
    then reads in the calculated solubility from "out.txt" and inserts this value into the
    database.
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
        ff.write('2 \n')       #number of components in system
        if (self['indices']['A']=='CO2'):
            ff.write('co2 \n')
        elif (self['indices']['A']=='Air'):
            ff.write('air \n')
        elif (self['indices']['A']=='CyP'):
            ff.write('cyclopentane \n')

        #pass molar liquid composition
        x1l_str = str(self['point']['xl1'])
        x2l_str = str(self['point']['xl2'])
        ff.write('%s \n' %(x1l_str))
        ff.write('%s \n' %(x2l_str))

        ff.close()

        #create output file for detailed code
        fff = open('out.txt', 'w+')
        fff.close()

        # Execute the detailed model
        # path to **this** file + /src/...
        # will break if distributed computing
        os.system(os.path.dirname(os.path.abspath(__file__))+\
            '/src/pcsaft')

        # Analyse output
        # os.getcwd() returns the path to the "launcher" directory
        try:
            FILE = open(os.getcwd()+'/out.txt','r')
        except IOError:
            raise IOError("File not found")
        self['point']['H'] = float(FILE.readline())
        FILE.close()

        return FWAction(mod_spec=[{'_push': self['point']}])


f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void surroSolubility
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

    const double term1 = P1*(1/T - 1/P2);
    const double term2 = exp(term1);

    outputs[0] = P0*term2;

    //outputs[0] = P0 + T*P1 + P2*T*T;
}
''',
    # These are global bounds for the function
    inputs={
        'T': { 'min': 200.0, 'max': 500.0},        #check if boundaries reasonable, from this range, the random values for the DOE are chosen!
        'xl1': { 'min': 0.0, 'max': 1.0 },
        'xl2': { 'min': 0.0, 'max': 1.0 },
    },
    outputs={
        'H': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -1E10, 'max': 1E10, 'argPos': 0 },    #check if boundaries are reasonable!!!
        'param1': { 'min': -1E10, 'max': 1E10, 'argPos': 1 },
        'param2': { 'min': -1.0e9, 'max': 1E10, 'argPos': 2 },
    },
    indices={
        'A': species
    },
)

m_solubilityCO2 = BackwardMappingModel(
    _id='Solubility[A=CO2]',
    surrogateFunction=f,
    exactTask=SolubilityExactSim(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 320, 350, 380],
            'xl1': [1.1e-3, 1.0e-3, 1.0e-3, 1.0e-4],
            'xl2': [0.9989, 0.999, 0.999, 0.9999],
        },
    ),
    outOfBoundsStrategy=Strategy.ExtendSpaceStochasticSampling(
        nNewPoints=4
    ),
    parameterFittingStrategy=Strategy.NonLinFitWithErrorContol(
        testDataPercentage=0.2,
        maxError=0.05,
        improveErrorStrategy=Strategy.StochasticSampling(
            nNewPoints=2
        ),
        maxIterations=5  # Currently not used
    ),
)
m_solubilityAir = BackwardMappingModel(
    _id='Solubility[A=Air]',
    surrogateFunction=f,
    exactTask=SolubilityExactSim(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 320, 350, 380],
            'xl1': [1.1e-3, 1.0e-3, 1.0e-3, 1.0e-4],
            'xl2': [0.9989, 0.999, 0.999, 0.9999],
        },
    ),
    outOfBoundsStrategy=Strategy.ExtendSpaceStochasticSampling(
        nNewPoints=4
    ),
    parameterFittingStrategy=Strategy.NonLinFitWithErrorContol(
        testDataPercentage=0.2,
        maxError=0.05,
        improveErrorStrategy=Strategy.StochasticSampling(
            nNewPoints=2
        ),
        maxIterations=5  # Currently not used
    ),
)
m_solubilityCyclopentane = BackwardMappingModel(
    _id='Solubility[A=CyP]',
    surrogateFunction=f,
    exactTask=SolubilityExactSim(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 320, 350, 380],
            'xl1': [1.1e-3, 1.0e-3, 1.0e-3, 1.0e-4],
            'xl2': [0.9989, 0.999, 0.999, 0.9999],
        },
    ),
    outOfBoundsStrategy=Strategy.ExtendSpaceStochasticSampling(
        nNewPoints=4
    ),
    parameterFittingStrategy=Strategy.NonLinFitWithErrorContol(
        testDataPercentage=0.2,
        maxError=0.05,
        improveErrorStrategy=Strategy.StochasticSampling(
            nNewPoints=2
        ),
        maxIterations=5  # Currently not used
    ),
)
