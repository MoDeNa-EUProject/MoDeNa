import os
import modena
from modena import ForwardMappingModel, BackwardMappingModel, SurrogateModel, \
CFunction
import modena.Strategy as Strategy
from fireworks.user_objects.firetasks.script_task import FireTaskBase, ScriptTask
from fireworks import Firework, Workflow, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from blessings import Terminal
from jinja2 import Template

## Create terminal for colour output
term = Terminal()




# For the case, when only foam conductivity and no aging is needed.
m = Strategy.BackwardMappingScriptTask(
        script=os.path.dirname(os.path.abspath(__file__))+'/src_dummy/workflowdummy'
            )
