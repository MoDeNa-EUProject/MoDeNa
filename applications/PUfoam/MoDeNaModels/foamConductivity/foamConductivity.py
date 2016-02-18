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
Surrogate function, model definition and backward mapping FireTask for
Foam conductivity model.

@author    Pavel Ferkl
@copyright 2014-2015, MoDeNa Project. GNU Public License.
@ingroup   app_aging
"""

import os
from modena import *
import modena.Strategy as Strategy
from fireworks.utilities.fw_utilities import explicit_serialize
from jinja2 import Template
import polymerConductivity
import gasConductivity
import gasMixtureConductivity


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
        f = open('foamConductivity.in', 'w')
        f.write('{0:.6e}\n'.format(temp+1))
        f.write('{0:.6e}\n'.format(temp-1))
        f.write('{0:.6e}\t{1:.6e}\t{2:.6e}\n'.format(xCO2,xAir,xCyP))
        f.write('0.9\n')
        f.write('0.9\n')
        f.write('1.2\n')
        f.write('1.1e3\n')
        f.write('{0:.6e}\n'.format(eps))
        f.write('{0:.6e}\n'.format(dcell))
        f.write('2\n')
        f.write('0.5e-6\n')
        f.write('{0:.6e}\n'.format(fstrut))
        f.write('1e-6\n')
        f.write('3e-2\n')
        f.write('200\n')
        f.write('10000\n')
        f.write('t\n')
        f.write('0.2\n')
        f.write('10\n')
        f.write('f\n')
        f.write('PeriodicRVEBoxStruts.vtk\n')
        f.close()
        # Execute the detailed model
        # path to **this** file + /src/...
        # will break if distributed computing
        os.system(os.path.dirname(os.path.abspath(__file__))+'/src/kfoam')
        # Analyse output
        # os.getcwd() returns the path to the "launcher" directory
        try:
            FILE = open(os.getcwd()+'/foamConductivity.out','r')
        except IOError:
            raise IOError("File not found")

        self['point']['kfoam'] = float(FILE.readline())

        os.remove('foamConductivity.in')
        os.remove('foamConductivity.out')

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

    fs=alpha*fstrut;
    Xs=(1+4*kgas/(kgas+polymer_thermal_conductivity))/3.0;
    Xw=2*(1+kgas/(2*polymer_thermal_conductivity))/3.0;
    X=(1-fs)*Xw+fs*Xs;
    kappa=4.09*sqrt(1-eps)/dcell;
    kr=16*sigma*pow(T,3)/(3*kappa);
    kfoam = (kgas*eps+polymer_thermal_conductivity*X*(1-eps))/(eps+(1-eps)*X)+beta*kr;

    outputs[0] = kfoam;
}
''',
    # These are global bounds for the function
    inputs={
        'eps': {'min': 0, 'max': 1},
        'dcell': {'min': 0, 'max': 1e-1},
        'fstrut': {'min': 0, 'max': 1},
        'gasMixtureConductivity': {'min': 0, 'max': 1e-1},
        'polymer_thermal_conductivity': {'min': 0, 'max': 1e0},
        'T': {'min': 273, 'max': 450},
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

# use input file to Foam aging application to initialize with reasonable data.
fname='foamAging.in'
try:
    f = open(os.getcwd()+'/../'+fname,'r')
except IOError:
    try:
        f = open(os.getcwd()+'/'+fname,'r')
    except IOError:
        f = open(os.getcwd()+'/example_inputs/'+fname,'r')

a=f.readline()
a=f.readline()
a=f.readline()
a=f.readline()
a=f.readline()
a=f.readline()
T0=float(a.split()[0])
a=f.readline()
rhop=float(a.split()[0])
a=f.readline()
a=f.readline()
a=f.readline()
xAir0=float(a.split()[0])
xCO20=float(a.split()[1])
xCyP0=float(a.split()[2])
a=f.readline()
a=f.readline()
a=f.readline()
dcell0=float(a.split()[0])
a=f.readline()
fstrut0=float(a.split()[0])
a=f.readline()
rho0=float(a.split()[0])
f.close()
eps0=1-rho0/rhop

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
