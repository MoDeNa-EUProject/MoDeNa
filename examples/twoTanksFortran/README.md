TWO TANKS EXAMPLE:
==================

The discharge of air from one tank into another through a nozzle. The
problem only makes sense, if you think of the flow through the nozzle as a
much more complex problem which you have to solve with - let's say 3D CFD. So
the two tanks are the macroscopic and the nozzle is the microscopic problem.
This is also a backward mapping problem if you assume that the range of
inputs is unknown a-priori. twoTanksMacroscopicProblem uses the MoDeNa
interface library to embed an even simpler model for the flow rate. While
twoTanksFullProblem implements it fully integrated.

TODO:
The example specific sources should be moved into the example, but this requires
building (and finding) the execuables here.


How to run?
-----------

# Make sure `PYTHONPATH` and `LD_LIBRARY_PATH` are set
# TODO:
# Make this easier to use
    export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${HOME}/lib/pkgconfig:/usr/local/lib/pkgconfig
    export PYTHONPATH=${PYTHONPATH:-}:${HOME}/lib/python2.7/site-packages
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib/python2.7/site-packages:${HOME}/lib:/usr/local/lib

# Compile project specific sources, i.e. "models":
    flowRate="../models/flowRate/src"
    twoTank="../models/twoTank/src"
    cmake -H${flowRate} -B${flowRate} && make --directory=${flowRate}
    cmake -H${twoTank} -B${twoTank} && make --directory=${twoTank}


# Compile project specific sources, i.e. "models":
    flowRate="../models/flowRate/src"
    twoTank="../models/twoTankFortran/src"
    cmake -H${flowRate} -B${flowRate} && make --directory=${flowRate}
    cmake -H${twoTank} -B${twoTank} && make --directory${twoTank}

# Initialise the model in the database
    ./initModel

# Start the workflow
    ./workflow

# Run again to see that no fitting is done on the second start
    ./workflow

