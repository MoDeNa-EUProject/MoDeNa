MULTICOMPONENT DIFFUSION EXAMPLE:
=================================

Binary diffusion by Fuller etal. is combined with Wilke multicomponent diffusion.

This example only shows the use of index sets. No backward mapping!


How to run?
-----------

# Make sure `PYTHONPATH` and `LD_LIBRARY_PATH` are set
# TODO:
# Make this easier to use
    export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${HOME}/lib/pkgconfig:/usr/local/lib/pkgconfig
    export PYTHONPATH=${PYTHONPATH:-}:${HOME}/lib/python2.7/site-packages
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib/python2.7/site-packages:${HOME}/lib:/usr/local/lib

# Compile project specific sources, i.e. "models":
    fuller="../models/fullerEtAlDiffusion/src"
    cmake -H${fuller} -B${fuller} && make --directory=${fuller}

# Initialise the model in the database
    ./initModel

# Start the workflow
    ./workflow

# Run again to see that no fitting is done on the second start
    ./workflow

