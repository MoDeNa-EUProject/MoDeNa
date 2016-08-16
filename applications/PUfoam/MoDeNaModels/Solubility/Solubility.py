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
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

import os
from modena import *
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
system = IndexSet(
    name= 'solubility_num_of_components',
    names= [ '2', '3' ]
)

@explicit_serialize
class SolubilityExactSim(ModenaFireTask):
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
    def task(self, fw_spec):
        # Write input for detailed model
        ff = open('in.txt', 'w')
        Tstr = str(self['point']['T'])
        ff.write('%s \n' %(Tstr))
        if (self['indices']['B']=='2'):
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
        elif (self['indices']['B']=='3'):
            ff.write('3 \n')       #number of components in system
            if (self['indices']['A']=='CO2'):
                ff.write('co2 \n')
            elif (self['indices']['A']=='Air'):
                ff.write('air \n')
            elif (self['indices']['A']=='CyP'):
                ff.write('cyclopentane \n')
            #pass molar liquid composition
            x1l_str = str(self['point']['xl1'])
            x2l_str = str(self['point']['xl2'])
            x3l_str = str(self['point']['xl3'])
            ff.write('%s \n' %(x1l_str))
            ff.write('%s \n' %(x2l_str))
            ff.write('%s \n' %(x3l_str))
        ff.close()

        #create output file for detailed code
        fff = open('out.txt', 'w+')
        fff.close()

        # Execute the detailed model
        # path to **this** file + /src/...
        # will break if distributed computing
        ret = os.system(os.path.dirname(os.path.abspath(__file__))+\
            '/src/pcsaft')
        # This call enables backward mapping capabilities
        self.handleReturnCode(ret)
        # Analyse output
        # os.getcwd() returns the path to the "launcher" directory
        try:
            FILE = open(os.getcwd()+'/out.txt','r')
        except IOError:
            raise IOError("File not found")
        self['point']['H'] = float(FILE.readline())
        FILE.close()


Ccode2='''
#include "modena.h"
#include "math.h"

void surroSolubility2
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

outputs[0] = P0 + T*P1 + P2*T*T;
}
'''
Ccode3='''
#include "modena.h"
#include "math.h"

void surroSolubility3
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

outputs[0] = P0 + T*P1 + P2*T*T;
}
'''
outputs={
    'H': { 'min': 9e99, 'max': -9e99, 'argPos': 0 }
}
parameters={
    'param0': { 'min': 1e-12, 'max': 1E1, 'argPos': 0 },    #check if boundaries are reasonable!!!
    'param1': { 'min': -1e-1, 'max': 1e-1, 'argPos': 1 },
    'param2': { 'min': -1e-1, 'max': 1e-1, 'argPos': 2 },
}
indices={
    'A': species,
    'B': system
}
inputs2={
    'T': { 'min': 200.0, 'max': 550.0},        #check if boundaries reasonable, from this range, the random values for the DOE are chosen!
    'xl1': { 'min': 0.0, 'max': 1.0 },
    'xl2': { 'min': 0.0, 'max': 1.0 },
}
inputs3={
    'T': { 'min': 200.0, 'max': 550.0},        #check if boundaries reasonable, from this range, the random values for the DOE are chosen!
    'xl1': { 'min': 0.0, 'max': 1.0 },
    'xl2': { 'min': 0.0, 'max': 1.0 },
    'xl3': { 'min': 0.0, 'max': 1.0 },
}
f2 = CFunction(Ccode=Ccode2,
    inputs=inputs2,
    outputs=outputs,
    parameters=parameters,
    indices=indices
)
f3 = CFunction(Ccode=Ccode3,
    inputs=inputs3,
    outputs=outputs,
    parameters=parameters,
    indices=indices
)

outOfBoundsStrategy=Strategy.ExtendSpaceStochasticSampling(
    nNewPoints=4
)
parameterFittingStrategy=Strategy.NonLinFitWithErrorContol(
    testDataPercentage=0.2,
    maxError=0.05,
    improveErrorStrategy=Strategy.StochasticSampling(
        nNewPoints=2
    ),
    maxIterations=5  # Currently not used
)
m_solubilityCO2PU = BackwardMappingModel(
    _id='Solubility[A=CO2,B=2]',
    surrogateFunction=f2,
    exactTask=SolubilityExactSim(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 320, 350, 380],
            'xl1': [1.1e-3, 1.0e-3, 1.0e-3, 1.0e-4],
            'xl2': [0.9989, 0.999, 0.999, 0.9999],
        },
    ),
    outOfBoundsStrategy=outOfBoundsStrategy,
    parameterFittingStrategy=parameterFittingStrategy
)
m_solubilityAirPU = BackwardMappingModel(
    _id='Solubility[A=Air,B=2]',
    surrogateFunction=f2,
    exactTask=SolubilityExactSim(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 320, 350, 380],
            'xl1': [1.1e-3, 1.0e-3, 1.0e-3, 1.0e-4],
            'xl2': [0.9989, 0.999, 0.999, 0.9999],
        },
    ),
    outOfBoundsStrategy=outOfBoundsStrategy,
    parameterFittingStrategy=parameterFittingStrategy
)
m_solubilityCyclopentanePU = BackwardMappingModel(
    _id='Solubility[A=CyP,B=2]',
    surrogateFunction=f2,
    exactTask=SolubilityExactSim(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 320, 350, 380],
            'xl1': [1.1e-3, 1.0e-3, 1.0e-3, 1.0e-4],
            'xl2': [0.9989, 0.999, 0.999, 0.9999],
        },
    ),
    outOfBoundsStrategy=outOfBoundsStrategy,
    parameterFittingStrategy=parameterFittingStrategy
)
# m_solubilityCO2 = BackwardMappingModel(
#     _id='Solubility[A=CO2,B=3]',
#     surrogateFunction=f3,
#     exactTask=SolubilityExactSim(),
#     substituteModels=[],
#     initialisationStrategy=Strategy.InitialPoints(
#         initialPoints={
#             'T': [290, 320, 350, 380],
#             'xl1': [0.43, 0.42, 0.40, 0.38],
#             'xl2': [0.53, 0.52, 0.50, 0.48],
#             'xl3': [0.04, 0.06, 0.10, 0.14],
#         },
#     ),
#     outOfBoundsStrategy=outOfBoundsStrategy,
#     parameterFittingStrategy=parameterFittingStrategy
# )
# m_solubilityAir = BackwardMappingModel(
#     _id='Solubility[A=Air,B=3]',
#     surrogateFunction=f3,
#     exactTask=SolubilityExactSim(),
#     substituteModels=[],
#     initialisationStrategy=Strategy.InitialPoints(
#         initialPoints={
#             'T': [280, 300, 320, 350],
#             'xl1': [0.08, 0.078, 0.074, 0.068],
#             'xl2': [0.85, 0.82, 0.79, 0.72],
#             'xl3': [0.07, 0.102, 0.136, 0.212],
#         },
#     ),
#     outOfBoundsStrategy=outOfBoundsStrategy,
#     parameterFittingStrategy=parameterFittingStrategy
# )
# m_solubilityCyclopentane = BackwardMappingModel(
#     _id='Solubility[A=CyP,B=3]',
#     surrogateFunction=f3,
#     exactTask=SolubilityExactSim(),
#     substituteModels=[],
#     initialisationStrategy=Strategy.InitialPoints(
#         initialPoints={
#             'T': [280, 300, 320, 350],
#             'xl1': [0.021, 0.023, 0.025, 0.027],
#             'xl2': [0.12e-4, 0.71e-4, 0.31e-3, 0.19e-2],
#             'xl3': [0.97886, 0.97629, 0.9719, 0.954],
#         },
#     ),
#     outOfBoundsStrategy=outOfBoundsStrategy,
#     parameterFittingStrategy=parameterFittingStrategy
# )
# below are experimental solubility models
# they are needed, because we want substitute model for bubble growth model
inputsExp={
    'T': { 'min': 200.0, 'max': 550.0}
}
CcodeR11Baser='''
#include "modena.h"
#include "math.h"

void surroSolubilityR11Baser
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
const double P3 = parameters[3];

outputs[0] = (P0 + P1*exp(-pow(T-P2,2)/(2*P3*P3)))*1100.0/137.37e-3/101e3;
}
'''
parameters4={
    'param0': { 'min': 1e-12, 'max': 1E1, 'argPos': 0 },
    'param1': { 'min': -1e-1, 'max': 1e-1, 'argPos': 1 },
    'param2': { 'min': -1e-1, 'max': 1e-1, 'argPos': 2 },
    'param3': { 'min': -1e-1, 'max': 1e-1, 'argPos': 3 },
}
fR11Baser = CFunction(Ccode=CcodeR11Baser,
    inputs=inputsExp,
    outputs=outputs,
    parameters=parameters4
)
## [R11, Baser](http://dx.doi.org/10.1002/pen.760340804)
parR11Baser = [1e-7, 4.2934, 203.3556, 40.016]
m_solubilityR11Baser = ForwardMappingModel(
    _id='SolubilityR11Baser',
    surrogateFunction=fR11Baser,
    substituteModels=[],
    parameters=parR11Baser,
    inputs=inputsExp,
    outputs=outputs,
)
CcodeCO2Baser='''
#include "modena.h"
#include "math.h"

void surroSolubilityCO2Baser
(
const modena_model_t* model,
const double* inputs,
double *outputs
)
{
{% block variables %}{% endblock %}

const double P0 = parameters[0];

outputs[0] = P0;
}
'''
parameters1={
    'param0': { 'min': 1e-12, 'max': 1E1, 'argPos': 0 },
}
if 'argPos' in inputsExp['T']: #added by previous CFunction
    inputsExp['T'].pop('argPos')
fCO2Baser = CFunction(Ccode=CcodeCO2Baser,
    inputs=inputsExp,
    outputs=outputs,
    parameters=parameters1
)
## [CO2, Baser](http://dx.doi.org/10.1002/pen.760340804)
parCO2Baser = [1.1e-4]
m_solubilityCO2Baser = ForwardMappingModel(
    _id='SolubilityCO2Baser',
    surrogateFunction=fCO2Baser,
    substituteModels=[],
    parameters=parCO2Baser,
    inputs=inputsExp,
    outputs=outputs,
)
CcodePentGupta='''
#include "modena.h"
#include "math.h"

void surroSolubilityPentGupta
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
const double P3 = parameters[3];
const double P4 = parameters[4];

outputs[0] = P0/(exp((P1-P2*T)/(P3-T))-P4)*1100.0/72.15e-3/101e3;
}
'''
parameters5={
    'param0': { 'min': 1e-12, 'max': 1E1, 'argPos': 0 },
    'param1': { 'min': -1e-1, 'max': 1e-1, 'argPos': 1 },
    'param2': { 'min': -1e-1, 'max': 1e-1, 'argPos': 2 },
    'param3': { 'min': -1e-1, 'max': 1e-1, 'argPos': 3 },
    'param4': { 'min': -1e-1, 'max': 1e-1, 'argPos': 4 },
}
if 'argPos' in inputsExp['T']: #added by previous CFunction
    inputsExp['T'].pop('argPos')
fPentGupta = CFunction(Ccode=CcodePentGupta,
    inputs=inputsExp,
    outputs=outputs,
    parameters=parameters5
)
## [n-pentane, Gupta](http://dx.doi.org/10.1002/pen.11405)
parPentGupta = [-3.3e-4, 2.09e4, 67.5, 8.69e4, 1.01]
m_solubilityPentGupta = ForwardMappingModel(
    _id='SolubilityPentGupta',
    surrogateFunction=fPentGupta,
    substituteModels=[],
    parameters=parPentGupta,
    inputs=inputsExp,
    outputs=outputs,
)
CcodePentWinkler='''
#include "modena.h"
#include "math.h"

void surroSolubilityPentWinkler
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
const double P3 = parameters[3];

outputs[0] = (P0 + P1*exp(-pow(T-P2,2)/(2*P3*P3)))*1100.0/72.15e-3/101e3;
}
'''
if 'argPos' in inputsExp['T']: #added by previous CFunction
    inputsExp['T'].pop('argPos')
fPentWinkler = CFunction(Ccode=CcodePentWinkler,
    inputs=inputsExp,
    outputs=outputs,
    parameters=parameters5
)
## n-pentane, Winkler
parPentWinkler = [0.0064,0.0551,298.0,17.8]
m_solubilityPentWinkler = ForwardMappingModel(
    _id='SolubilityPentWinkler',
    surrogateFunction=fPentWinkler,
    substituteModels=[],
    parameters=parPentWinkler,
    inputs=inputsExp,
    outputs=outputs,
)

@explicit_serialize
class SolubilityExactSim2(ModenaFireTask):
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
    def task(self, fw_spec):
        # Write input for detailed model
        ff = open('in.txt', 'w')
        Tstr = str(self['point']['T'])
        ff.write('%s \n' %(Tstr))
        if (self['indices']['B']=='2'):
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
        elif (self['indices']['B']=='3'):
            ff.write('3 \n')       #number of components in system
            if (self['indices']['A']=='CO2'):
                ff.write('co2 \n')
            elif (self['indices']['A']=='Air'):
                ff.write('air \n')
            elif (self['indices']['A']=='CyP'):
                ff.write('cyclopentane \n')
            #pass molar liquid composition
            # x1l_str = str(self['point']['xl1'])
            # x2l_str = str(self['point']['xl2'])
            # x3l_str = str(self['point']['xl3'])
            x1l_str = 1.0e-4
            x2l_str = 0.5
            x3l_str = 0.5
            # x1l_str = 0.2
            # x2l_str = 0.4
            # x3l_str = 0.4
            ff.write('%s \n' %(x1l_str))
            ff.write('%s \n' %(x2l_str))
            ff.write('%s \n' %(x3l_str))
        ff.close()

        #create output file for detailed code
        fff = open('out.txt', 'w+')
        fff.close()

        # Execute the detailed model
        # path to **this** file + /src/...
        # will break if distributed computing
        ret = os.system(os.path.dirname(os.path.abspath(__file__))+\
            '/src/pcsaft')
        # This call enables backward mapping capabilities
        self.handleReturnCode(ret)

        # Analyse output
        # os.getcwd() returns the path to the "launcher" directory
        try:
            FILE = open(os.getcwd()+'/out.txt','r')
        except IOError:
            raise IOError("File not found")
        self['point']['H'] = float(FILE.readline())
        FILE.close()

Ccode='''
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
const double P3 = parameters[3];

const double term1 = P1*(1/T - 1/P2);
const double term2 = exp(term1);

outputs[0] = P0*term2;

outputs[0] = P0 + T*P1 + P2*T*T;

outputs[0] = P0+P1*exp(-pow((T-P2),2)/P3);
}
'''
if 'argPos' in inputsExp['T']: #added by previous CFunction
    inputsExp['T'].pop('argPos')
parameters4={
    'param0': { 'min': -9e9, 'max': 9e9+2*1.42885440e-04, 'argPos': 0 },    #check if boundaries are reasonable!!!
    'param1': { 'min': -9e9, 'max': 9e9+2*1.22132172e+00, 'argPos': 1 },
    'param2': { 'min': -9e9, 'max': 9e9+2*-7.97789449e+02, 'argPos': 2 },
    'param3': { 'min': -9e9, 'max': 9e9+2*1.33835999e+05, 'argPos': 3 },
}
fx = CFunction(Ccode=Ccode,
    inputs=inputsExp,
    outputs=outputs,
    parameters=parameters4,
    indices=indices
)
outOfBoundsStrategy=Strategy.ExtendSpaceStochasticSampling(
    nNewPoints=4
)
parameterFittingStrategy=Strategy.NonLinFitWithErrorContol(
    testDataPercentage=0.2,
    maxError=0.05,
    improveErrorStrategy=Strategy.StochasticSampling(
        nNewPoints=2
    ),
    maxIterations=5  # Currently not used
)
m_solubilityCO2 = BackwardMappingModel(
    _id='Solubility[A=CO2,B=3]',
    surrogateFunction=fx,
    exactTask=SolubilityExactSim2(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [290, 450, 350, 380, 550],
        },
    ),
    outOfBoundsStrategy=outOfBoundsStrategy,
    parameterFittingStrategy=parameterFittingStrategy
)
m_solubilityCyclopentane = BackwardMappingModel(
    _id='Solubility[A=CyP,B=3]',
    surrogateFunction=fx,
    exactTask=SolubilityExactSim2(),
    substituteModels=[],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints={
            'T': [300, 350, 380, 420, 400, 450, 280, 550],
        },
    ),
    outOfBoundsStrategy=outOfBoundsStrategy,
    parameterFittingStrategy=parameterFittingStrategy
)
