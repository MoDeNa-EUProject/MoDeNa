'''

   ooo        ooooo           oooooooooo.             ooooo      ooo
   `88.       .888'           `888'   `Y8b            `888b.     `8'
    888b     d'888   .ooooo.   888      888  .ooooo.   8 `88b.    8   .oooo.
    8 Y88. .P  888  d88' `88b  888      888 d88' `88b  8   `88b.  8  `P  )88b
    8  `888'   888  888   888  888      888 888ooo888  8     `88b.8   .oP"888
    8    Y     888  888   888  888     d88' 888    .o  8       `888  d8(  888
   o8o        o888o `Y8bod8P' o888bood8P'   `Y8bod8P' o8o        `8  `Y888""8o

Copyright
    2014 MoDeNa Consortium, All rights reserved.

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
    Christos Mitrias
'''

import os
import modena
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, CFunction
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

# ********************************* Class ************************************ #
@explicit_serialize
class rheologyExactTask(FireTaskBase):
    """
    A FireTask that starts a microscopic code and updates the database.
    """   
    def run_task(self, fw_spec):
        print(
            term.yellow +
            "Performing exact simulation (microscopic code recipe)" +
            term.normal
        )


        # Write input
        f = open('rheologyExact.in', 'w')

        f.close()

        # Execute the detailed model
        os.system('../src/rheologyexact')

        # Analyse output
        f=open('rheologyExact.out','r')
        self['point']['mu'] = float(f.readline())
        f.close()
        
        os.remove('rheologyExact.in')
        os.remove('rheologyExact.out')

        return FWAction(mod_spec=[{'_push': self['point']}])
       
f = CFunction(
    Ccode= '''
#include "modena.h"
#include "math.h"

void rheology_SM
(
    const double* parameters,
    const double* inherited_inputs,
    const double* inputs,
    double *outputs
)
{    
    const double temp = inputs[0]; // temperature
    const double shear = inputs[1]; // shear rate
    const double conv = inputs[2]; // conversion
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
    const double conv_gel = 0.615;
    const double mu_0_const = 0.195; 
    const double mu_inf_const = 0.266;
//    const double lambda = 11.35 ;
//    const double alpha = 2; 
//    const double n_rh = 0.2;

    double mu_0, mu_inf, f_t;
    double mu_ap;
    
    mu_0 = (log(conv+mu_d) - log(mu_d) + pow(conv_gel / ( conv_gel - conv ), mu_a + conv*mu_b + mu_c*pow(conv,2))) * mu_0_const; 
    mu_inf = (log(conv+mu_d) - log(mu_d) + pow(conv_gel / ( conv_gel - conv ), mu_a + conv*mu_b + mu_c*pow(conv,2)))* mu_inf_const;
    f_t = A_mu * exp(E_mu / R_rh / temp );

    mu_ap = (mu_inf + (mu_0 - mu_inf)*pow(1 + pow(lambda*shear,alpha), (n_rh - 1) / alpha)) * f_t;

    outputs[0] = mu_ap;
}
''',
   # These are global bounds for the function
   inputs={
       'temp': {'min': 0, 'max': 9e99, 'argPos': 0 },
       'shear': {'min': 0, 'max': 9e99, 'argPos': 1 },
       'conv': {'min': 0, 'max': 1, 'argPos': 2 },
   },
   outputs={
       'mu_ap': { 'min': 0, 'max': 9e99, 'argPos': 0 },
   },
   parameters={
       'lamdba': { 'min': 11.35, 'max': 11.35, 'argPos': 0 },
       'alpha': { 'min': 2, 'max': 2, 'argPos': 1 },
       'n_rh': { 'min': 0.2, 'max': 0.2, 'argPos': 2 },
   },
)

m = ForwardMappingModel(
    _id= 'rheology',
    surrogateFunction= f,
    substituteModels= [ ],
    parameters= [ 11.35, 2, 0.2],
    inputs={
        'temp': { 'min': 0, 'max': 9e99 },
        'shear': { 'min': 0, 'max': 9e99 },
        'conv': { 'min': 0, 'max': 9e99 },
    },
    outputs={
        'mu_ap': {'min': 0, 'max': 9e99 },
    },
)

