#!/usr/bin/python
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
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.
@endcond'''
from __future__ import division
"""
@file
Creates input files for Bubble growth and QmomKinetics detailed models.

@author Pavel Ferkl
@copyright 2014-2015, MoDeNa Project. GNU Public License.
"""

import json
from numpy import exp,log,pi
# read json
f=open("unifiedInput.json",'r')
inputs=json.load(f)
f.close()
# prepare input for Bubble growth model
f=open("inputs.in",'w')
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["integrator"]))
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["method"]))
f.write('\n')
if (inputs["bubbleGrowth"]["inertialTerm"]):
    f.write('t\n')
else:
    f.write('f\n')
if (inputs["bubbleGrowth"]["solubilityCorrection"]):
    f.write('t\n')
else:
    f.write('f\n')
f.write('{0:.3g}\n'.format(inputs["bubbleGrowth"]["meshCoarseningParameter"]))
f.write('\n')
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["internalNodes"]))
f.write('{0:.3g}\n'.format(inputs["bubbleGrowth"]["initialTime"]))
f.write('{0:.3g}\n'.format(inputs["bubbleGrowth"]["finalTime"]))
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["outerTimeSteps"]))
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["maxInnerTimeSteps"]))
f.write('{0:.3g}\n'.format(inputs["bubbleGrowth"]["relativeTolerance"]))
f.write('{0:.3g}\n'.format(inputs["bubbleGrowth"]["absoluteTolerance"]))
f.write('\n')
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["numberOfDissolvedGases"]))
f.write('{0:d}\n'.format(inputs["bubbleGrowth"]["carbonDioxidePosition"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["pressure"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["molarMass"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["molarMass"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["heatCapacity"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["heatCapacityInLiquidPhase"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["heatCapacityInLiquidPhase"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["heatCapacityInGaseousPhase"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["heatCapacityInGaseousPhase"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["evaporationHeat"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["evaporationHeat"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["PBL"]["density"]))
f.write('\n')
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["temperature"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["bubbleRadius"]))
NN=inputs["initialConditions"]["numberBubbleDensity"]
R0=inputs["initialConditions"]["bubbleRadius"]
f.write('{0:.3g}\n'.format((1/(NN*exp(log(4/3*pi*R0**3)))+1)**(1/3)))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["concentrations"]["polyol"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["concentrations"]["water"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["concentrations"]["isocyanate"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["initialConditions"]["concentrations"]["blowingAgents"]["PBL"],\
    inputs["initialConditions"]["concentrations"]["blowingAgents"]["CO2"]))
f.write('1 0 0\n')
f.write('\n')
if (inputs["kinetics"]["kineticModel"]=="Baser"):
    f.write('1\n')
elif (inputs["kinetics"]["kineticModel"]=="modena"):
    f.write('2\n')
elif (inputs["kinetics"]["kineticModel"]=="BaserRx"):
    f.write('3\n')
else:
    print 'kinetic model unknown in Bubble growth model'
    exit
if (inputs["kinetics"]["useDilution"]):
    f.write('t\n')
else:
    f.write('f\n')
f.write('{0:.3g}\n'.format(inputs["kinetics"]["gellingReaction"]["frequentialFactor"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["gellingReaction"]["activationEnergy"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["blowingReaction"]["frequentialFactor"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["blowingReaction"]["activationEnergy"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["gellingReaction"]["reactionEnthalpy"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["blowingReaction"]["reactionEnthalpy"]))
f.write('\n')
f.write('{0:d}\n'.format(inputs["physicalProperties"]["polymer"]["polymerDensityModel"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["density"]))
f.write('\n')
f.write('{0:d}\n'.format(inputs["physicalProperties"]["surfaceTensionModel"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["surfaceTension"]))
f.write('\n')
f.write('{0:d} {1:d}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["diffusivityModel"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["diffusivityModel"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["diffusivity"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["diffusivity"]))
f.write('\n')
f.write('{0:d} {1:d}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["solubilityModel"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["solubilityModel"]))
f.write('{0:.3g} {1:3g}\n'.format(\
    inputs["physicalProperties"]["blowingAgents"]["PBL"]["solubility"],\
    inputs["physicalProperties"]["blowingAgents"]["CO2"]["solubility"]))
f.write('\n')
f.write('{0:d}\n'.format(inputs["physicalProperties"]["polymer"]["viscosityModel"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["viscosity"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["maxViscosity"]))
f.write('4.1e-8\n')
f.write('38.3e3\n')
f.write('0.85\n')
f.write('4e0\n')
f.write('-2e0\n')
f.close()
# prepare input for QmomKinetics
f=open("inputsQmom.in",'w')
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["pressure"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["temperature"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["gellingReaction"]["frequentialFactor"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["gellingReaction"]["activationEnergy"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["concentrations"]["polyol"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["concentrations"]["isocyanate"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["concentrations"]["water"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["blowingReaction"]["frequentialFactor"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["blowingReaction"]["activationEnergy"]))
f.write('8.3145\n')
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["density"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["PBL"]["density"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["gellingReaction"]["reactionEnthalpy"]))
f.write('{0:.3g}\n'.format(inputs["kinetics"]["blowingReaction"]["reactionEnthalpy"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["heatCapacity"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["CO2"]["heatCapacityInLiquidPhase"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["PBL"]["heatCapacityInGaseousPhase"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["PBL"]["heatCapacityInLiquidPhase"]))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["PBL"]["evaporationHeat"]))
if (inputs["physicalBlowingAgent"]=="pentane"):
    f.write('1\n')
elif (inputs["physicalBlowingAgent"]=="R11"):
    f.write('2\n')
else:
    print 'unknown blowing agent (solubility)'
    exit
f.write('2\n') #TODO implement
if (inputs["kinetics"]["kineticModel"]=="Baser"):
    f.write('1\n')
elif (inputs["kinetics"]["kineticModel"]=="BaserRx"):
    f.write('2\n')
else:
    print 'kinetic model unknown in QmomKinetics'
    exit
if (inputs["kinetics"]["useDilution"]):
    f.write('1\n')
else:
    f.write('0\n')
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["CO2"]["molarMass"]*1e3))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["blowingAgents"]["PBL"]["molarMass"]*1e3))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["polymer"]["molarMassNCO"]*1e3))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["air"]["molarMass"]*1e3))
H=inputs["physicalProperties"]["blowingAgents"]["CO2"]["solubility"]
p=inputs["physicalProperties"]["pressure"]
M=inputs["physicalProperties"]["blowingAgents"]["CO2"]["molarMass"]
rho=inputs["physicalProperties"]["polymer"]["density"]
f.write('{0:.3g}\n'.format(H*p*M/rho))
c=inputs["initialConditions"]["concentrations"]["blowingAgents"]["PBL"]
rho=inputs["physicalProperties"]["polymer"]["density"]
M=inputs["physicalProperties"]["blowingAgents"]["PBL"]["molarMass"]
f.write('{0:.3g}\n'.format(c*M/rho))
c=inputs["initialConditions"]["concentrations"]["blowingAgents"]["CO2"]
rho=inputs["physicalProperties"]["polymer"]["density"]
M=inputs["physicalProperties"]["blowingAgents"]["CO2"]["molarMass"]
f.write('{0:.3g}\n'.format(c*M/rho))
f.write('{0:.3g}\n'.format(inputs["physicalProperties"]["surfaceTension"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["bubbleRadiusDeviation"]))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["bubbleRadius"]*2))
f.write('{0:.3g}\n'.format(inputs["initialConditions"]["numberBubbleDensity"]))
f.write('{0:.3g}\n'.format(inputs["QmomKinetics"]["absoluteTolerance"]))
f.write('{0:.3g}\n'.format(inputs["QmomKinetics"]["relativeTolerance"]))
f.write('{0:.3g}\n'.format(inputs["QmomKinetics"]["timeStep"]))
f.write('{0:.3g}\n'.format(inputs["QmomKinetics"]["endTime"]))
