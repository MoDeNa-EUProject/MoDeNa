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
'''

import os
import modena
from modena import BackwardMappingModel
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal

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
class ModenaBackwardMappingTask(ScriptTask):
    """
    A FireTask that starts a macroscopic code and catches its return code.
    @author: Henrik Rusche
    """
    required_params = ['script']

    def run_task(self, fw_spec):
        print(
            term.yellow
          + "Performing backward mapping simulation (macroscopic code recipe)"
          + term.normal
        )

        self['defuse_bad_rc'] = True

        # Execute the macroscopic code by calling function in base class
        ret = super(ModenaBackwardMappingTask, self).run_task(fw_spec)

        # Analyse return code

        print('return code = %i' % ret.stored_data['returncode'])
        if ret.stored_data['returncode'] > 199:
            print term.cyan + "Performing Design of Experiments" + term.normal
            # TODO
            # Finding the 'failing' model using the outsidePoint will fail
            # eventually fail when running in parallel. Need to pass id of
            # calling FireTask. However, this requires additional code in the
            # library as well as cooperation of the recipie
            model = modena.SurrogateModel.loadFailing()

            return model.outOfBoundsFwAction(
                model,
                self,
                outsidePoint= model.outsidePoint
            )

        else:
            print('We are done')
            return ret

m = ModenaBackwardMappingTask(
    script= os.path.dirname(os.path.realpath(__file__)) + \
            '/src/twoTanksMacroscopicProblem'
)
