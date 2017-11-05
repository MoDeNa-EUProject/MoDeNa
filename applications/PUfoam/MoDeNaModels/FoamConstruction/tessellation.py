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
Uses Neper 3.
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
    centers_rads = np.loadtxt(myfile12, usecols=(0, 1, 2, 3))
    centers_rads[:, 3] = centers_rads[:, 3] / 2.0
    centers = centers_rads[:, :3]  # All centers of spheres
    rads = centers_rads[:, 3]  # All radii of spheres
    max_centers = max(centers.tolist())
    min_centers = min(centers.tolist())
    edge_cube_size = [math.ceil(max_centers[0] - min_centers[0]),
                      math.ceil(max_centers[1] - min_centers[1]),
                      math.ceil(max_centers[2] - min_centers[2])]
    edge_rve_size = int(max(edge_cube_size))  # For NEPER: Size of edge of RVE
    # create periodic RVE directly using Neper's new periodicity option
    myfile3 = os.path.join(MYPATH, 'centers.txt')
    with open(myfile3, 'w') as fout:
        for i in range(0, 1):
            for j in range(0, number_of_cells):
                fout.write('{0:f}\t{1:f}\t{2:f}\n'.format(centers[j][i],
                                                          centers[j][i + 1],
                                                          centers[j][i + 2]))
    myfile4 = os.path.join(MYPATH, 'Rads.txt')
    with open(myfile4, 'w') as fout:
        for i in range(0, 1):
            for j in range(0, number_of_cells):
                fout.write('{0:f}\n'.format(rads[j]))
    # Note: Neper regularization is not available for periodic tessellations
    command_tessellation = "neper -T \
        -n {0:d} \
        -domain 'cube({1:d},{2:d},{3:d})' \
        -periodicity x,y,z \
        -morpho voronoi \
        -morphooptiini 'coo:file(centers.txt),weight:file(Rads.txt)' \
        -o {4} -format tess,geo \
        -statcell vol -statedge length -statface area \
        -statver x".format(number_of_cells, edge_rve_size, edge_rve_size,
                           edge_rve_size, filename)
    os.system(command_tessellation)
    if visualize_tessellation:  # needs POV-Ray
        command_visualization = "neper -V {0}RVE27.tess -datacellcol ori \
            -datacelltrs 0.5 -showseed all -dataseedrad @Rads.txt \
            -dataseedtrs 1.0 -print {0}RVE27".format(filename)
        os.system(command_visualization)
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
