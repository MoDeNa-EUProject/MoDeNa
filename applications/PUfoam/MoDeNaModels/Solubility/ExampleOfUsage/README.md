
Example simulation with Solubility model:

A simple macroscopic code calls the solubility model for different temperatures.


How to run?
-----------

# Make sure PYTHONPATH and LD_LIBRARY_PATH are set
# TODO:
# Make this easier to use
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${HOME}/lib/pkgconfig:/usr/local/lib/pkgconfig
export PYTHONPATH=${PYTHONPATH:-}:${HOME}/lib/python2.7/site-packages
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib/python2.7/site-packages:${HOME}/lib:/usr/local/lib

# Compile project specific sources (only locally)
cd src
cmake .
make


# Initialise the model in the database
./initModel

# Start the workflow
./workflow
