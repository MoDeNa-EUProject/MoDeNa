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
@file       SurfaceTension.py
@namespace  SurfaceTension.SurfaceTension
@ingroup    mod_surfacetension
@brief      Surrogate Function, Surrogate Model templates and Model Recipe
@author    Jonas Mairhofer
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@details

# Surface Tension python module.

Contains the following:

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

"""

import os
from jinja2 import Template
from fireworks.utilities.fw_utilities import explicit_serialize
from modena import CFunction, BackwardMappingModel, IndexSet, ModenaFireTask
from modena import Strategy

## @var blowing_agents
# @brief (MoDeNa) Index Set for the Blowing Agents that are valid for the model
# @details
#
# The index set contains two elements:
#
# @f[
#     \mathbb{A} = \left\{ \text{Air}, \text{CO2} \right\}
# @f]
#
blowing_agents = IndexSet(
    name= 'blowing_agents',
    names= [ 'AIR', 'CO2']
)

## @var monomers
# @brief (MoDeNa) Index Set for the Monomers that are valid for the model
# @details
#
# The index set contains four elements:
#
# @f[
#   \mathbb{B} = \left\{ \text{PU}, \text{THF}, \text{HEXANE} \right\}
# @f]
#
monomers = IndexSet(
    name = 'monomers',
    names = ['PU', 'THF', 'HEXANE']
)

## @var surfactant
# @brief (MoDeNa) Index Set for the Surfactants that are valid for the model
# @details
#
# The index set contains four elements:
#
# @f[
#   \mathbb{C}= \left\{ \text{surfactant}, \text{no\_surfactant} \right\}
# @f]
#
surfactant = IndexSet(
    name = 'surfactant',
    names = ['surfactant','no_surfactant']
)



# ********************************* Class ********************************** #
@explicit_serialize
class SurfaceTensionExactSim(ModenaFireTask):
    """
    @brief    Recipe for PCSAFT Surface Tension Application
    @details
              This FireTask controls the execution of the detailed model of the
              Surface Tension model. The detailed model is a density functional
              theory implementation based on PC-SAFT. A detailed description of
              this model can be found in Deliverable 1.3 on the [MoDeNa]
              website.

              In order to start the detailed model, the input values for the
              model are first written to the file "in.txt". The detailed model
              code picks them up from this file and performs the according
              calculation. Once it is done, the output value is written to the
              file "out.txt". This FireTask then reads in the calculated
              surface tension from "out.txt" and inserts this value into the
              database.

              [MoDeNa]: http://www.modenaproject.eu/

    @attention
    @pre
    @note
    @warning
    """

    def task(self, fw_spec):

        # Write input for detailed model
        self.generate_inputfile()

        # Execute detailed model
        run_command = os.path.dirname(os.path.abspath(__file__))+'/src/PCSAFT_SurfaceTension -snes_monitor_short -ksp_monitor_short \
                   -nx 800 -rc 9.0 -box 300 -erel 1e-08 -init_pert 0 \
                   -snes_type newtonls  -snes_converged_reason  \
                   -snes_atol 1e-07 -snes_rtol 1e-07 -snes_stol 1e-07 -snes_max_it 20 \
                   -ksp_max_it 15 -ksp_gmres_restart 50 \
                   -snes_linesearch_type l2 -snes_linesearch_damping 0.3 -snes_linesearch_monitor \
                   -snes_max_fail 1 -snes_max_linear_solve_fail 100 \
                   -ksp_gmres_cgs_refinement_type refine_always \
                   -snes_ksp_ew -snes_ksp_ew_version 1 -snes_ksp_ew_rtol0 0.5 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_threshold 0.1 \
                   -jac 0 -pc_type none > log'


        ret = os.system(run_command)
        # This call enables backward mapping capabilities (not needed in this example)
        self.handleReturnCode(ret)

        # Analyse output
        self.analyse_output()

    def generate_inputfile(self):
        """
        @brief   Generate a input file.
        """
        with open('in.txt','w') as f:
            f.write("{}\n".format(self['point']['T']))
            if self['indices']['C']=='no_surfactant':
                ncomp=2
                f.write("{}\n".format(ncomp))
                f.write("{}\n".format(self['indices']['A'].lower()))
                f.write("{}\n".format(self['indices']['B'].lower()))
            else:
                ncomp=3
                f.write("{}\n".format(ncomp))
                f.write("{}\n".format(self['indices']['A'].lower()))
                f.write("{}\n".format(self['indices']['B'].lower()))
                f.write("{}\n".format(self['indices']['C'].lower()))
            for i in range(ncomp):
                f.write("{}\n".format(0.0))
        with open('out.txt','w+') as FILE:
            pass

    def analyse_output(self):
        """ analysing the output of the file.
            @TODO consider adding check for empty file
        """
        with open('out.txt', 'r') as FILE:
            self['point']['ST'] = float(FILE.readline())

## @var f
# @brief (MoDeNa) Surrogate Function Template
# @details
#
# The surrogate function is an indexed function of the form:
#
# @f[
#    \hat{\mathcal{M}} := f_{A,B,C}(T; \theta_1, \theta_2, \theta_3) \quad;\quad 270 \leq T \leq 550 \quad A\in\mathbb{A} \quad B\in \mathbb{B} \quad C\in \mathbb{C}
# @f]
#
# Where @f$\mathbb{A}@f$ and @f$\mathbb{B}@f$ respectively are the index sets @ref blowing_agents and @ref monomers.
#
# @f[
#    f_{A,B}(T; \theta_1, \theta_2, \theta_3) := \theta_1 \cdot T + \theta_2 \cdot T^2 + \theta_3 \cdot T^3
# @f]
#
# @attention Look here
# @pre
# @note
# @todo     Check if parameter bounds are reasonable
# @warning
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
    inputs={
        'T': { 'min': 270.0, 'max': 550.0 },
    },
    outputs={
        'ST': { 'min': 9e99, 'max': -9e99, 'argPos': 0 },
    },
    parameters={
        'param0': { 'min': -1E10, 'max': 1E10, 'argPos': 0 },
        'param1': { 'min': -1E10, 'max': 1E10, 'argPos': 1 },
        'param2': { 'min': -1E10, 'max': 1E10, 'argPos': 2 },
    },
    species = {
        'A' : blowing_agents,
        'B' : monomers,
        'C' : surfactant,
    }
)

## @var m
# @brief (MoDeNa) Surrogate Model Template (Air, THF)
# @details
#
# The surrogate model is defined for the indices "Air" and "THF":
# @f[
#    \hat{\mathcal{M}} := f_{\text{Air},\text{THF},\text{surfactant}}( T; \theta_1, \theta_2, \theta_3)
# @f]
#
# @attention Look here
# @pre
# @note
# @todo
# @warning
# @bug
#
m = BackwardMappingModel(
    _id= 'SurfaceTension[A=AIR,B=PU,C=surfactant]',
    surrogateFunction= f,
    exactTask= SurfaceTensionExactSim(),
    substituteModels= [ ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [290.0, 300.0, 315.0, 350.0, 375.0, 400.0, 450.0, 500.0],
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

## @var m2
# @brief (MoDeNa) Surrogate Model Template (Air, PU)
# @details
#
# @f[
#    \hat{\mathcal{M}} := f_{\text{Air},\text{PU},\text{no\_surfactant}}( T; \theta_1, \theta_2, \theta_3)
# @f]
#
# @note
# @todo
# @warning
# @bug
#
m2 = BackwardMappingModel(
    _id= 'SurfaceTension[A=AIR,B=PU,C=no_surfactant]',
    surrogateFunction= f,
    exactTask= SurfaceTensionExactSim(),
    substituteModels= [ ],
    initialisationStrategy= Strategy.InitialPoints(
        initialPoints=
        {
            'T': [290.0, 300.0, 315.0, 350.0, 375.0, 400.0, 450.0, 500.0],
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
