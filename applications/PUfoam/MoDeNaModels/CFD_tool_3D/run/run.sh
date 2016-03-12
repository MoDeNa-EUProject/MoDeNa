#!/bin/bash
# export PKG_CONFIG_PATH=${PKG_CONFIG_PATH:-}:${HOME}/lib/pkgconfig:/usr/local/lib/pkgconfig
# export PYTHONPATH=${PYTHONPATH:-}:${HOME}/lib/python2.7/site-packages/modena
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:${HOME}/lib/python2.7/site-packages:${HOME}/lib/modena:/usr/local/lib
cp -r ../0 ./
cp -r ../constant ./
cp -r ../system ./
cp ../getMeReady.sh ./
cp ../case.foam ./
./getMeReady.sh
blockMesh
setFields
rm -fv log
echo "PUFoam is running..."
${FOAM_USER_APPBIN}/PUFoam >& log
echo "Done."