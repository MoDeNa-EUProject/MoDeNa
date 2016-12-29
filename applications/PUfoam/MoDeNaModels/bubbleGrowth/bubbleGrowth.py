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
@file       bubbleGrowth.py
@namespace  bubbleGrowth.bubbleGrowth
@ingroup    mod_bubbleGrowth
@brief      Backward mapping Firetask for Bubble growth model.
@author     Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@details

# Bubble growth python module.

Contains a FireTask, which runs the detailed model and copies the results into
the results folder. The relative path and name of the detailed model executable
are hard coded. The FireTask is BackwardMapping, which means that if one of the
lower scale backward mapping models will get out of validity range that model
will be re-fitted to larger range and the detailed model will be re-run. This is
repeated until the detailed model succesfully finishes (or possibly crashes for
other reason, in which case an error is printed out).
"""

import os
from modena.Strategy import BackwardMappingScriptTask

## @var m
# @brief Bubble Growth Application Recipe
# @details
# Runs the detailed model and saves the results.
m = BackwardMappingScriptTask(
    script=os.path.dirname(os.path.abspath(__file__))
    + '/src/bblgrExact'
    + ' && cp *.out *.txt ../results/bubbleGrowth/'
)
