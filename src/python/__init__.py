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

    The Modena interface library is free software; you can redistribute it
    and/or modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    Modena is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with Modena.  If not, see <http://www.gnu.org/licenses/>.
@endcond'''

"""
@file
Module providing the MoDeNa python interface

@copyright  2014-2016, MoDeNa Project. GNU Public License.
"""

import os, sys
from pkg_resources import get_distribution

__version__ = get_distribution('modena').version

MODENA_INSTALL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODENA_WORKING_DIR = os.path.realpath(os.getcwd())

from Strategy import BackwardMappingScriptTask, ModenaFireTask
from SurrogateModel import CFunction, IndexSet, \
    SurrogateModel, ForwardMappingModel, BackwardMappingModel, \
    ModenaFireTask, MODENA_PARSED_URI

def find_module(target, startsearch=MODENA_WORKING_DIR):
    """Function recursively searching through the file tree for "target"

    @arg target: 'str' name of directory e.g. "Desktop"
    """

    pth = os.path.abspath(startsearch)
    while target not in os.listdir(pth):
        pth = os.path.abspath(os.path.join(pth,'..'))  # step back a directory
        if os.path.ismount(pth):                       # break if we hit "root"
            pth = None
            break

    if pth is not None:
        sys.path.insert(0, os.path.join(pth,target))
    else:
        print "Could not find directory: %s" %(target)

def import_helper():
    from os.path import dirname
    import imp


    fp = None
    try:
        fp, pathname, description = imp.find_module(
            'libmodena',
            [ find_module("modena", os.path.join(MODENA_INSTALL_DIR, "..","..")) ]
        )
    except ImportError:
        import libmodena
        return libmodena
    if fp is not None:
        try:
            _mod = imp.load_module('libmodena', fp, pathname, description)
        finally:
            fp.close()
        return _mod


find_module("MoDeNaModels")   # Look for a models directory
libmodena = import_helper()
del import_helper

##
# @defgroup python_interface_library
# Module providing the MoDeNa python interface

