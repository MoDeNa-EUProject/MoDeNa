"""
@file       skeleton.py
@namespace  FoamConstruction.skeleton
@ingroup    mod_foamConstruction
@brief      Packing and tesselation.
@author     Mohammad Marvi-Mashhadi
@author     Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@details
Prepares representative volume element (RVE) of foam.
"""
import os
import time
import random
import math
import numpy as np
import geo_tools
# current directory
MYPATH = os.getcwd()


def pack_spheres(average_radius, deviance, number_of_cells, alternative_algorithm):
    """Packs spheres into a periodic domain. Creates Project01.rco with sphere
    centers and radii."""
    Rad = abs((np.random.normal(average_radius, deviance, number_of_cells)))
    Rads1 = list(range(number_of_cells))
    t = 0
    for i in range(number_of_cells):
        c = abs(Rad[t])
        Rads1[t] = c.astype(np.float)
        t = t + 1
    Rads1 = sorted(Rads1)
    v = 0.00
    for i in range(number_of_cells):
        v = v + ((2.00 * Rads1[i])**3.00)
    centers = [[0 for i in range(3)] for j in range(number_of_cells)]
    v = v * 1.40
    lc = v**(1.00 / 3.00)
    K = 0
    while K == 0:
        j = -1
        h = 0
        timeout = time.time() + 10
        while number_of_cells >= j and h == 0:
            if time.time() > timeout:
                h = 1
                break
            j = j + 1
            if j == number_of_cells:
                K = 1
                break
            PickCenterX, PickCenterY, PickCenterZ =\
                lc * random.random(),\
                lc * random.random(),\
                lc * random.random()
            while (lc - Rads1[j] >= PickCenterX
                    and lc - Rads1[j] >= PickCenterY
                    and lc - Rads1[j] >= PickCenterZ and Rads1[j] < PickCenterX
                    and Rads1[j] < PickCenterY and Rads1[j] < PickCenterZ):
                PickCenterX, PickCenterY, PickCenterZ =\
                    lc * random.random(),\
                    lc * random.random(),\
                    lc * random.random()
            centers[j][0], centers[j][1], centers[j][2] =\
                PickCenterX, PickCenterY, PickCenterZ
            KeepCentreX, KeepCentreY, KeepCentreZ, KeepR =\
                PickCenterX, PickCenterY, PickCenterZ, Rads1[j]
            if j > 0:
                for t in range(0, j):
                    if ((((((KeepCentreX - centers[t][0])**2.00) +
                            ((KeepCentreY - centers[t][1])**2.00) +
                            ((KeepCentreZ - centers[t][2])**2.00))**0.50) -
                            (KeepR + Rads1[t])) < 0.000) and t != j:
                        centers[j][0], centers[j][0], centers[j][0] = 0, 0, 0
                        j = j - 1
                        break
    myfile = os.path.join(MYPATH, 'Project01.rco')
    f = open(myfile, 'w')
    for i in range(number_of_cells):
        f.write('{0:f}{1}{2:f}{3}{4:f}{5}{6:f}\n'.format(
            centers[i][0], '	', centers[i][1], '	', centers[i][2], '	',
            2.0 * Rads1[i]))
    f.close()
    if alternative_algorithm:
        os.system('./spherepack')


def tessellate(filename, number_of_cells, visualize_tessellation):
    """Use Laguerre tessellation from Neper to create dry foam. Uses Project01.rco
    as input file."""
    myfile12 = os.path.join(MYPATH, 'Project01.rco')
    CentersRads = np.loadtxt(myfile12, usecols=(0, 1, 2, 3))
    Centers = CentersRads[:, :3]  # All centers of spheres
    Rads = CentersRads[:, 3]  # All radii of spheres
    Rads = Rads / 2
    MAXcenters = max(Centers.tolist())
    Mincenters = min(Centers.tolist())
    EdgeCubeSize = [math.ceil(MAXcenters[0] - Mincenters[0]),
                    math.ceil(MAXcenters[1] - Mincenters[1]),
                    math.ceil(MAXcenters[2] - Mincenters[2])]
    EdgeRVESize = int(max(EdgeCubeSize))  # For NEPER: Size of edge of RVE
    EdgeRVESize27 = 3 * EdgeRVESize
    MaxCenters = np.amax(Centers, axis=0)
    L2 = int(0.5 + MaxCenters[0])
    X = L2 + Centers[:, 0]
    Y = L2 + Centers[:, 1]
    Z = L2 + Centers[:, 2]
    # translate seeds in 26 directions to simulate periodicity
    Centers27 = np.array([
        X, X + L2, X - L2, X + L2, X, X + L2, X, X - L2, X - L2,
        X - L2, X + L2, X, X, X - L2, X + L2, X + L2, X, X,
        X - L2, X, X, X - L2, X + L2, X + L2, X - L2, X - L2, X + L2,
        Y, Y + L2, Y - L2, Y + L2, Y + L2, Y, Y - L2, Y, Y - L2,
        Y + L2, Y - L2, Y - L2, Y + L2, Y, Y, Y, Y + L2, Y,
        Y, Y - L2, Y, Y + L2, Y - L2, Y + L2, Y - L2, Y + L2, Y - L2,
        Z, Z + L2, Z - L2, Z, Z + L2, Z + L2, Z - L2, Z - L2, Z,
        Z, Z, Z + L2, Z - L2, Z + L2, Z - L2, Z, Z, Z + L2,
        Z, Z, Z - L2, Z + L2, Z + L2, Z - L2, Z + L2, Z - L2, Z - L2
    ])
    # Creation of input .txt files for neper
    myfile3 = os.path.join(MYPATH, 'Centers27.txt')
    ff = open(myfile3, 'w')
    for i in range(0, 27):
        for j in range(0, number_of_cells):
            ff.write('{0:f}\t{1:f}\t{2:f}\n'.format(Centers27[i, j],
                                                    Centers27[i + 27, j],
                                                    Centers27[i + 54, j]))
    ff.close()
    myfile4 = os.path.join(MYPATH, 'Rads27.txt')
    fff = open(myfile4, 'w')
    for i in range(0, 27):
        for j in range(0, number_of_cells):
            fff.write('{0:f}\n'.format(Rads[j]))
    fff.close()
    commandTessellation = "neper -T \
        -n {0:d} \
        -domain 'cube({1:d},{2:d},{3:d})' \
        -morpho voronoi \
        -morphooptiini 'coo:file(Centers27.txt),weight:file(Rads27.txt)' \
        -o {4}RVE27 -format tess,geo \
        -statcell vol -statedge length -statface area \
        -statver x".format(27 * number_of_cells, EdgeRVESize27, EdgeRVESize27,
                           EdgeRVESize27, filename)
    os.system(commandTessellation)
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
        -o {4}Closed -format tess,geo \
        -statcell vol -statedge length -statface area \
        -statver x".format(number_of_cells, EdgeRVESize, EdgeRVESize,
                           EdgeRVESize, filename)
    os.system(commandTessellation)
    if visualize_tessellation:  # needs POV-Ray
        commandVisualization = "neper -V {0}RVE27.tess -datacellcol ori \
            -datacelltrs 0.5 -showseed all -dataseedrad @Rads.txt \
            -dataseedtrs 1.0 -print {0}RVE27".format(filename)
        os.system(commandVisualization)
    sdat = geo_tools.read_geo(filename + ".geo")
    edat = geo_tools.extract_data(sdat)
    point = edat["point"]
    line = edat["line"]
    with open('{0}.gnu'.format(filename), 'w') as flp:
        for pidx in line.itervalues():
            flp.write('{0} {1} {2}\n'.format(
                point[pidx[0]][0], point[pidx[0]][1], point[pidx[0]][2]))
            flp.write('{0} {1} {2}\n\n\n'.format(
                point[pidx[1]][0], point[pidx[1]][1], point[pidx[1]][2]))
