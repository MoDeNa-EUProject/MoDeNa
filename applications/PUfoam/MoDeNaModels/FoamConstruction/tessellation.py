"""
@file       tessellation.py
@namespace  FoamConstruction.tessellation
@ingroup    mod_foamConstruction
@brief      Tesselation.
@author     Mohammad Marvi-Mashhadi
@author     Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@details
Prepares representative volume element (RVE) of foam using Laguerre tessellation.
"""
import os
import math
import numpy as np
from geo_tools import read_geo, extract_data
# current directory
MYPATH = os.getcwd()


def tessellate(filename, number_of_cells, visualize_tessellation):
    """Use Laguerre tessellation from Neper to create dry foam. Uses Project01.rco
    as input file."""
    myfile12 = os.path.join(MYPATH, 'Project01.rco')
    CentersRads = np.loadtxt(myfile12, usecols=(0, 1, 2, 3))
    CentersRads[:, 3] = CentersRads[:, 3] / 2.0
    Centers = CentersRads[:, :3]  # All centers of spheres
    Rads = CentersRads[:, 3]  # All radii of spheres
    Rads = Rads / 2
    MAXcenters = max(Centers.tolist())
    Mincenters = min(Centers.tolist())
    EdgeCubeSize = [math.ceil(MAXcenters[0] - Mincenters[0]),
                    math.ceil(MAXcenters[1] - Mincenters[1]),
                    math.ceil(MAXcenters[2] - Mincenters[2])]
    EdgeRVESize = int(max(EdgeCubeSize))  # For NEPER: Size of edge of RVE
    # create periodic RVE directly using Neper's new periodicity option
    myfile3 = os.path.join(MYPATH, 'Centers.txt')
    ff = open(myfile3, 'w')
    for i in range(0, 1):
        for j in range(0, number_of_cells):
            ff.write('{0:f}\t{1:f}\t{2:f}\n'.format(Centers[j][i],
                                                    Centers[j][i + 1],
                                                    Centers[j][i + 2]))
    ff.close()
    myfile4 = os.path.join(MYPATH, 'Rads.txt')
    fff = open(myfile4, 'w')
    for i in range(0, 1):
        for j in range(0, number_of_cells):
            fff.write('{0:f}\n'.format(Rads[j]))
    fff.close()
    commandTessellation = "neper -T \
        -n {0:d} \
        -domain 'cube({1:d},{2:d},{3:d})' \
        -periodicity x,y,z \
        -morpho voronoi \
        -morphooptiini 'coo:file(Centers.txt),weight:file(Rads.txt)' \
        -o {4} -format tess,geo \
        -statcell vol -statedge length -statface area \
        -statver x".format(number_of_cells, EdgeRVESize, EdgeRVESize,
                           EdgeRVESize, filename)
    os.system(commandTessellation)
    if visualize_tessellation:  # needs POV-Ray
        commandVisualization = "neper -V {0}RVE27.tess -datacellcol ori \
            -datacelltrs 0.5 -showseed all -dataseedrad @Rads.txt \
            -dataseedtrs 1.0 -print {0}RVE27".format(filename)
        os.system(commandVisualization)
    sdat = read_geo(filename + ".geo")
    edat = extract_data(sdat)
    point = edat["point"]
    line = edat["line"]
    with open('{0}.gnu'.format(filename), 'w') as flp:
        for pidx in line.itervalues():
            flp.write('{0} {1} {2}\n'.format(
                point[pidx[0]][0], point[pidx[0]][1], point[pidx[0]][2]))
            flp.write('{0} {1} {2}\n\n\n'.format(
                point[pidx[1]][0], point[pidx[1]][1], point[pidx[1]][2]))
