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
Surrogate function, model definition and backward mapping FireTask for
Foam conductivity model.

@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
from modena import *
import modena.Strategy as Strategy
from fireworks.utilities.fw_utilities import explicit_serialize
from jinja2 import Template
import json
import polymerConductivity
import gasConductivity
import gasMixtureConductivity
from json import encoder
encoder.FLOAT_REPR = lambda o: format(o, '.12g')

@explicit_serialize
class FoamConductivityExactTask(ModenaFireTask):
    """
    A FireTask that starts a microscopic code and updates the database.
    """
    def task(self, fw_spec):
        eps = self['point']['eps']
        dcell = self['point']['dcell']
        fstrut = self['point']['fstrut']
        temp = self['point']['T']
        xCO2 = self['point']['x[CO2]']
        xCyP = self['point']['x[CyP]']
        xO2 = self['point']['x[O2]']
        xN2 = self['point']['x[N2]']
        xAir = xN2+xO2
        # Write input
        inputs={"upperBoundary": {"temperature": temp+1,"emittance": 0.9}}
        inputs["lowerBoundary"]={"temperature": temp-1,"emittance": 0.9}
        inputs["gasComposition"]={"CO2": xCO2,"Air": xAir,"Cyclopentane": xCyP}
        inputs["gasDensity"]=1.2
        inputs["solidDensity"]=1.1e3
        inputs["porosity"]=eps
        inputs["cellSize"]=dcell
        inputs["morphologyInput"]=2
        inputs["wallThickness"]=0.5e-6
        inputs["strutContent"]=fstrut
        inputs["strutSize"]=1e-6
        inputs["foamThickness"]=3e-2
        inputs["spatialDiscretization"]=200
        inputs["useWallThicknessDistribution"]=True
        inputs["wallThicknessStandardDeviation"]=0.2
        inputs["numberOfGrayBoxes"]=10
        inputs["numericalEffectiveConductivity"]=False
        # inputs["structureName"]=
        inputs["testMode"]=False
        with open('foamConductivity.json','w') as f:
            json.dump(inputs, f, indent=4)
        # Execute the detailed model
        # path to **this** file + /src/...
        # will break if distributed computing
        ret = os.system(os.path.dirname(os.path.abspath(__file__))+'/src/kfoam')
        # This call enables backward mapping capabilities
        self.handleReturnCode(ret)
        # Analyse output
        # os.getcwd() returns the path to the "launcher" directory
        try:
            FILE = open(os.getcwd()+'/foamConductivity.out','r')
        except IOError:
            raise IOError("File not found")

        self['point']['kfoam'] = float(FILE.readline())

        # os.remove('foamConductivity.json')
        # os.remove('foamConductivity.out')

## Surrogate function for thermal conductivity of the foam.
#
# Foam conductivity is a function of porosity, cell size, strut content,
# conductivity of gas and solid phase and temperature.
f_foamConductivity = CFunction(
    Ccode='''
#include "modena.h"
#include "math.h"

void tcfoam_SM
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double alpha = parameters[1];
    const double beta = parameters[0];

    const double sigma=5.67e-8;

    double fs,Xs,Xw,X,kappa,kr;
    double kfoam;
    double kgas=gasMixtureConductivity;
    double kpol=polymer_thermal_conductivity;

    fs=alpha*fstrut;
    Xs=(1+4*kgas/(kgas+kpol))/3.0;
    Xw=2*(1+kgas/(2*kpol))/3.0;
    X=(1-fs)*Xw+fs*Xs;
    kappa=4.09*sqrt(1-eps)/dcell;
    kr=16*sigma*pow(T,3)/(3*kappa);
    kfoam = (kgas*eps+kpol*X*(1-eps))/(eps+(1-eps)*X)+beta*kr;

    outputs[0] = kfoam;
}
''',
    # These are global bounds for the function
    inputs={
        'eps': {'min': 0, 'max': 0.995},
        'dcell': {'min': 0, 'max': 1e-1},
        'fstrut': {'min': 0, 'max': 1},
        'gasMixtureConductivity': {'min': 0, 'max': 1e-1},
        'polymer_thermal_conductivity': {'min': 0, 'max': 1e0},
        'T': {'min': 273, 'max': 550},
        'x': {'index': gasConductivity.species, 'min': 0, 'max': 1},
    },
    outputs={
        'kfoam': {'min': 0, 'max': 1e0, 'argPos': 0},
    },
    parameters={
        'param1': {'min': -1e9, 'max': 1e9 + 2, 'argPos': 0},
        'param2': {'min': -1e9, 'max': 1e9 + 2, 'argPos': 1},
    },
)

# When initializing for Foam aging
# use input file to Foam aging application to initialize with reasonable data.
try:
    with open('foamAging.json','r') as f:
        inputs=json.load(f)
        T0=inputs['foamCondition']['conductivityTemperature']
        rhop=inputs['physicalProperties']['polymerDensity']
        xAir0=inputs['foamCondition']['initialComposition']['Air']
        xCO20=inputs['foamCondition']['initialComposition']['CO2']
        xCyP0=inputs['foamCondition']['initialComposition']['Cyclopentane']
        dcell0=inputs['morphology']['cellSize']
        fstrut0=inputs['morphology']['strutContent']
        rho0=inputs['morphology']['foamDensity']
        eps0=1-rho0/rhop
except IOError: # set some dummy values when not initializing
    eps0=0.96
    dcell0=300e-6
    fstrut0=0.8
    T0=283
    xCO20=0.0
    xCyP0=0.99
    xAir0=0.0
# set initial points so that they are close to one value
def setIP(a0):
    a=[]
    for i in xrange(4):
        a.append(a0)
    upar=1-1e-4
    opar=1+1e-4
    a[0]=a[0]*upar
    if a0==0:
        a[1]=1e-4
    else:
        a[1]=a[1]*opar
    return a

initialPoints_foamConductivity_auto = {
    'eps': setIP(eps0),
    'dcell': setIP(dcell0),
    'fstrut': setIP(fstrut0),
    'T': setIP(T0),
    'x[CO2]': setIP(xCO20),
    'x[CyP]': setIP(xCyP0),
    'x[O2]': setIP(xAir0*0.21),
    'x[N2]': setIP(xAir0*0.79),
}
# when testing,
# initialize for any composition to avoid getting out of bounds too many times
test=False
if test:
    initialPoints_foamConductivity_auto['x[CO2]']=[1,0,0,0]
    initialPoints_foamConductivity_auto['x[CyP]']=[0,1,0,0]
    initialPoints_foamConductivity_auto['x[O2]']=[0,0,1,0]
    initialPoints_foamConductivity_auto['x[N2]']=[0,0,0,1]
# When initializing for Foam expansion
# use dummy data to initialize to avoid getting out of bounds
foaming_ini={
    'eps': [0.9,0.0,0.96,0.99,0.7,0.5],
    'dcell': [200e-6,0.0,300e-6,100e-6,1e-2,200e-6],
    'fstrut': [0.0,1.0,0.7,0.6,0.0,0.9],
    'T': [280,549,300,350,330,300],
    'x[CO2]': [0,0,1,1,1,0],
    'x[CyP]': [1,1,0,0,0,1],
    'x[O2]': [0,0,0,0,0,0],
    'x[N2]': [0,0,0,0,0,0],
}
# we are initializing for Foam expansion if unifiedInput.json exists
try:
    with open('inputs/unifiedInput.json','r') as f:
        inputs=json.load(f)
    initialPoints_foamConductivity_auto = foaming_ini
except IOError:
    try:
        with open('../inputs/unifiedInput.json','r') as f:
            inputs=json.load(f)
        initialPoints_foamConductivity_auto = foaming_ini
    except:
        pass
## Surrogate model for foam conductivity
#
# Backward mapping model is used.
m_foamConductivity = BackwardMappingModel(
    _id='foamConductivity',
    surrogateFunction=f_foamConductivity,
    exactTask=FoamConductivityExactTask(),
    substituteModels=[
        gasMixtureConductivity.m_gasMixtureConductivity,\
        polymerConductivity.m_polymer_thermal_conductivity\
    ],
    initialisationStrategy=Strategy.InitialPoints(
        initialPoints=initialPoints_foamConductivity_auto,
    ),
    outOfBoundsStrategy=Strategy.ExtendSpaceStochasticSampling(
        nNewPoints=4
    ),
    parameterFittingStrategy=Strategy.NonLinFitWithErrorContol(
        testDataPercentage=0.2,
        maxError=0.01,
        improveErrorStrategy=Strategy.StochasticSampling(
            nNewPoints=2
        ),
        maxIterations=5  # Currently not used
    ),
)

## Foam conductivity simulation
#
# For the case, when only foam conductivity and no aging is needed.
m_simulation = Strategy.BackwardMappingScriptTask(
    script=os.path.dirname(os.path.abspath(__file__))+'/src/kfoam' +
        ' && cp foamConductivity.out ../results/' +
        ' && cp hahtf.out ../results/'
)
