#!/usr/bin/python
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
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.

Description
    A simple workflow

Authors
    Henrik Rusche

Contributors
'''

import modena
from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire
from modulefinder import ModuleFinder

from modena.Strategy import BackwardMappingScriptTask
import os

# Source code in src/twoTanksMacroscopicProblem.C
SIMULATION = BackwardMappingScriptTask(
    script='../run.sh'
)


# set up the LaunchPad and reset it
launchpad = LaunchPad()
launchpad.reset('', require_password=False)

# create the individual FireWorks and Workflow
# Source code in src/twoTanksMacroscopicProblem.C
wf = Workflow([Firework(SIMULATION)], {}, name="simulation")

# store workflow and launch it locally
launchpad.add_wf(wf)
rapidfire(launchpad)

