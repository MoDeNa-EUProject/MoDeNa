#!/bin/bash

# delete everything but the good stuff!
ls > list
egrep -v 'CMakeCache.txt|CMakeFiles|cmake_install.cmake|Makefile|CMakeLists.txt|cleanupScript.sh|coalescence.h|eigen|experimentalInputs.h|getBoost.sh|gnuplot_script.gnu|growth.h|initializeMoments.h|liquidBA.h|partialPressure.h|README|readParameters.h|results|testingKinetics|pda.h|QmomKinetics|QmomKinetics.cpp|write_kinetics.h|momentsConverter.h|plot_cfd.vsz|inputsQmom.in' list > list2
mv list2 list
rm -rf $(<list)
