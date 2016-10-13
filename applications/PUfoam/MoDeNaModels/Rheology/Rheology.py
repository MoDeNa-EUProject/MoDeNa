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
@author    Christos Mitrias
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   app_foaming
"""

from os import getcwd, makedirs, remove, system, symlink, walk
from os.path import abspath, dirname, isfile, join, relpath
from shutil import copy2
from glob import glob
import re
import json
from numpy.linalg import pinv
from numpy import sin, cos, pi, matrix, array, dot

import modena
import SurfaceTension
import polymerViscosity 
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, CFunction, ModenaFireTask
import modena.Strategy as Strategy

from fireworks.utilities.fw_utilities import explicit_serialize
from jinja2 import Template
import json

from math import pi,sqrt
import numpy as np

# ----------------------- Convenience variables ----------------------------- #
MODULE_DIR = dirname(abspath(__file__))
EXECUTION_DIR = getcwd()

with open(getcwd()+'/inputs/unifiedInput.json') as jsonfile:
    inputs=json.load(jsonfile)
    X_gel=inputs['kinetics']['gelPoint']

# ------------------------- Utility functions ------------------------------- #
def link_files(src, dst):
    """ Function symlinking all files and directories from directory "src"
        to the directory "dst".
    """
    for dirpath,_,filenames in walk(src):
        dstpth = join(dst, dirpath[len(src)+1:])
        relpth = relpath(dirpath, dstpth)
    
        #print "In ", dirpath
        #print "dstPth = ", dstpth
        #print "relpath = ", relpth
        #print filenames
    
        try:
            makedirs(dstpth)
        except:
            pass
    
        lines = []
        try:
            with open(join(dirpath,".filesToCopy")) as FILE:
                lines = FILE.read().splitlines()
        except:
            pass
    
        for fn in filenames:
            if any(fnmatch.fnmatch(fn, gl) for gl in lines):
                dstfn = join(dstpth, fn)
                srcfn = join(dirpath, fn)
                if not isfile(dstfn):
                    #print "copying", fn
                    copy2(srcfn, dstfn)
            else:
                dstfn = join(dstpth, fn)
                srcfn = join(relpth, fn)
                if not isfile(dstfn):
                    #print "linking", fn
                    symlink(srcfn, dstfn)


# -------------- Surrogate Functions AND Exact Tasks ------------------------ #

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Dummy model
@explicit_serialize
class RheologyExactTask_dummy(ModenaFireTask):
    """
    A FireTask that starts a microscopic code and updates the database.
    """
    INPUT_FILE = 'RheologyExact.in'
    OUTPUT_FILE = 'RheologyExact.out'
    APPLICATION_CMD = join(MODULE_DIR ,'src/rheologyexactdummy')

    def task(self, fw_spec):

        # Write input
        Template('''{{ s['point']['T'] }}
                {{ s['point']['shear'] }}
                {{ s['point']['X'] }}
                {{ s['point']['m0'] }}
                {{ s['point']['m1'] }}
		{{ s['point']['mu'] }}
		{{ s['point']['ST'] * 0,001 }}'''.strip()).stream(s=self).dump(self.INPUT_FILE)


        # Execute the detailed model
        ret = system(self.APPLICATION_CMD)
        # This enables backward mapping capabilities (not needed in this example)
        self.handleReturnCode(ret)

        # Analyse output
        with open(self.OUTPUT_FILE,'r') as FILE:
            self['point']['mu_car'] = float(FILE.readline())

        remove(self.INPUT_FILE)
        remove(self.OUTPUT_FILE)


f_dummy = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"
#include <stdio.h>

void rheology_dummy_SM
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{
    {% block variables %}{% endblock %}

    const double lambda = parameters[0];
    const double alpha = parameters[1];
    const double n_rh = parameters[2];

    const double A_mu = 0.0387; // could be loaded from some library for consistence
    const double E_mu = 10000;
    const double R_rh = 8.314;
    const double mu_a = 1.5;
    const double mu_b = 1;
    const double mu_c = 0;
    const double mu_d = 0.001;
    const double X_gel = '''+str(X_gel)+''';
    const double mu_0_const = 0.195;
    const double mu_inf_const = 0.266;

    double mu_0, mu_inf;
    double mu_car;

    if (X<X_gel) {
        mu_0 = (log(X+mu_d) - log(mu_d) + pow(X_gel / ( X_gel - X ), mu_a + X*mu_b + mu_c*pow(X,2))) * mu_0_const;
        mu_inf = (log(X+mu_d) - log(mu_d) + pow(X_gel / ( X_gel - X ), mu_a + X*mu_b + mu_c*pow(X,2)))* mu_inf_const;

        mu_car = (mu_inf + (mu_0 - mu_inf)*pow(1 + pow(lambda*shear,alpha), (n_rh - 1) / alpha));
    } else {
        mu_car = 1e6;
    }
    //    printf("apparent viscosity car %f", mu_car);
    outputs[0] = mu_car;
}
''',
   # These are global bounds for the function
   inputs={
       'T': {'min': 0, 'max': 550 },
       'shear': {'min': 0, 'max': 9e99 },
       'X': {'min': 0, 'max': 1 },
       'm0': {'min': 0, 'max': 9e99 },
       'm1': {'min': 0, 'max': 9e99 },
       'mu': {'min': 0, 'max': 1000 },
       'ST': {'min': 0, 'max': 100 },
   },
   outputs={
       'mu_car': { 'min': 0, 'max': 9e99, 'argPos': 0 },
   },
   parameters={
       'lamdba': { 'min': 1.35, 'max': 21.35, 'argPos': 0 },
       'alpha': { 'min': 0, 'max': 2, 'argPos': 1 },
       'n_rh': { 'min': 0, 'max': 2, 'argPos': 2 },
   },
)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> tFEM model
@explicit_serialize
class RheologyExactTask_tFEM(ModenaFireTask):
    """
    A FireTask that starts a microscopic code and updates the database.
    """
    inputTemplate = Template('''{{ s['point']['T'] }}
                {{ s['point']['shear'] }}
                {{ s['point']['X'] }}
                {{ s['point']['mu'] }}
                {{ s['point']['ST'] * 0.001 }}
                {{ s['point']['m0'] }}	
                {{ s['point']['m1'] }}'''.strip())

    def post_processing(self, fname="RheologyExact.out"):
        """Analyse output
        """
        fname = glob("bulk_stress_*.out")
        print "Reading output file: %s" %(fname)
        # 1)                                          Read and process raw data
        fc = tuple(open(fname[0], 'r'))#                              Read file
        DATA =  [[float(n) for n in re.findall(r'-?\d+.?\d+',l)] for l in fc]
        del fc#                                         Delete raw file content

        # 2)          First line of file contains parameters for the simulation
        header = ('strain', 'omega', 'dt', 'eta_p', 'Gamma', 'r')
        params = { k: v for (k, v) in zip(header, DATA.pop(0)) }
        # DATA structure: [[t, Gxx, Gxy, Gyy, dt, eta, dot_gamma, gamma]_1,...]

        # 3)                       Remove effects due to the initial conditions
        Tc = 2*pi/params['omega']/params['dt']
        tc = int(1*Tc)
        te = int(6*Tc)

        K = (te-tc)/1
        DATA = zip(*DATA[tc:te])

        # DATA = [ col for col in zip(*DATA) ]#                 Remove entries
        # DATA structure now transposed: [[t_0 ... t_n], [Gxx_1 ... Gxx_n],...]

        # 4)                               Generate input and response matrices
        Y = np.matrix( DATA[2] ).transpose()#                      Response matrix

        si = params['strain']*np.sin(params['omega']*np.array(DATA[0]) )#  Sine terms
        co = params['strain']*np.cos(params['omega']*np.array(DATA[0]) )# Cosine term
        X = np.matrix( zip(*[si, co]) )#                              Input matrix

        del si, co, DATA

        # 5)                                          Perform linear regression
        b = np.linalg.pinv(X)*Y#                       pinv: Moore-Penrose pseudo-inverse

        # 6)                                  eta_posity = abs(b0 + b1) / omega
        bsum = np.dot( b.T, b ).A1[0]#               Inner product sums "b" vector
        eta = sqrt( bsum )/params['omega']#

        return eta

    def execute_simulation(self):
        """ Method executing the detailed model
        """
        link_files(join(MODULE_DIR,"src"), getcwd())
        return system(MODULE_DIR+'/src/rheology')

    def task(self, fw_spec):
        """Method executed by FireWorks """

        fin="RheologyExact.in"

        # This enables backward mapping capabilities (not needed in this example)
        self.inputTemplate.stream(s=self).dump(fin)

        # call execution
        self.handleReturnCode( self.execute_simulation() )

        # call post processing
        self['point']['mu_car'] = self.post_processing()



f_tFEM = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"
#include <stdio.h>

void rheology_tFEM_SM
(
    const modena_model_t* model,
    const double* inputs,
    double *outputs
)
{    
    {% block variables %}{% endblock %}

    const double lambda = parameters[0];
    const double alpha = parameters[1];
    const double n_rh = parameters[2];

    const double A_mu = 0.0387; // could be loaded from some library for consistence
    const double E_mu = 10000;
    const double R_rh = 8.314;
    const double mu_a = 1.5;
    const double mu_b = 1;
    const double mu_c = 0;
    const double mu_d = 0.001;
    const double X_gel = 0.615;
    const double mu_0_const = 0.195; 
    const double mu_inf_const = 0.266;
//    const double lambda = 11.35 ;
//    const double alpha = 2;
//    const double n_rh = 0.2;

    double mu_0, mu_inf, f_t;
    double mu_car;

    mu_0 = (log(X+mu_d) - log(mu_d) + pow(X_gel / ( X_gel - X ), mu_a + X*mu_b + mu_c*pow(X,2))) * mu_0_const; 
    mu_inf = (log(X+mu_d) - log(mu_d) + pow(X_gel / ( X_gel - X ), mu_a + X*mu_b + mu_c*pow(X,2)))* mu_inf_const;
    f_t = A_mu * exp(E_mu / R_rh / T );

    mu_car = (mu_inf + (mu_0 - mu_inf)*pow(1 + pow(lambda*shear,alpha), (n_rh - 1) / alpha)) * f_t;
    //printf("apparent viscosity %f", mu_ap);
    outputs[0] = mu_car;
}
''',
   # These are global bounds for the function
   inputs={
       'T': {'min': 0, 'max': 550 },
       'shear': {'min': 0, 'max': 9e99 },
       'X': {'min': 0, 'max': 1 },
       'm0' : {'min': 0, 'max' : 9e99},
       'm1' : {'min': 0, 'max' : 9e99},
       'mu': {'min': 0, 'max': 1000 },
       'ST': {'min': 0, 'max': 100 },
   },
   outputs={
       'mu_car': { 'min': 0, 'max': 9e99, 'argPos': 0 },
   },
   parameters={
       'lamdba': { 'min': 1.35, 'max': 21.35, 'argPos': 0 },
       'alpha': { 'min': 0, 'max': 2, 'argPos': 1 },
       'n_rh': { 'min': 0, 'max': 2, 'argPos': 2 },
   },
)

# --------------------------------------------------------------------------- #

m = BackwardMappingModel(
    _id= 'Rheology',
    surrogateFunction= f_dummy,
    exactTask= RheologyExactTask_dummy(),
#    surrogateFunction= f_tFEM,
#    exactTask= RheologyExactTask_tFEM(),
    substituteModels= [ polymerViscosity.m_polymerViscosity, SurfaceTension.m],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [270.0, 330.0],
            'shear': [0, 10000],
            'X': [0, 0.3],
            'm0': [ 2.5e10, 1.1e10],
            'm1': [ 0.1, 0.1],
        },
    ),
    outOfBoundsStrategy= Strategy.ExtendSpaceStochasticSampling(
        nNewPoints= 4
    ),
    parameterFittingStrategy= Strategy.NonLinFitWithErrorContol(
        testDataPercentage= 0.2,
        maxError= 0.5,
        improveErrorStrategy= Strategy.StochasticSampling(
            nNewPoints= 2
        ),
        maxIterations= 5 # Currently not used
    ),
)

