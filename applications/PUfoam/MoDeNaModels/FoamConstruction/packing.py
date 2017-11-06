"""
@file       packing.py
@namespace  FoamConstruction.packing
@ingroup    mod_foamConstruction
@brief      Sphere packing.
@author     Mohammad Marvi-Mashhadi
@author     Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@details
Prepares packed spheres for tessellation.
"""
from __future__ import division, print_function
import struct
import os
import time
import random
import subprocess
import numpy as np
from numpy import array, pi, linspace, exp, log, sqrt
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import spack


def simple_packing(shape, scale, number_of_cells):
    "Simple and fast algorithm for packing"
    Rad = lognorm.rvs(shape, scale=scale, size=number_of_cells)
    print(Rad)
    Rad /= 2
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
    data = zip(*centers)
    data.append(Rads1)
    data = np.array(zip(*data))
    data[:, 3] = 2 * data[:, 3]
    return data


def create_input(npart, domain=1.0):
    """Create input file for packing program."""
    txt = """Particles count: {0}
Packing size: {1} {1} {1}
Generation start: 1
Seed: 341
Steps to write: 1000
Boundaries mode: 1
Contraction rate: 1.328910e-005
1. boundaries mode: 1 - bulk; 2 - ellipse (inscribed in XYZ box, Z is length of an ellipse); 3 - rectangle
2. generationMode = 1 (Poisson, R) or 2 (Poisson in cells, S)
    """.format(npart, domain)
    with open('generation.conf', 'w') as fout:
        fout.write(txt)


def make_csd(shape, scale, npart, show_plot=False):
    """Create cell size distribution and save it to file."""
    if shape == 0:
        rads = [scale + 0 * x for x in range(npart)]
    else:
        rads = lognorm.rvs(shape, scale=scale, size=npart)
    with open('diameters.txt', 'w') as fout:
        for rad in rads:
            fout.write('{0}\n'.format(rad))
    if shape == 0:
        xpos = linspace(scale / 2, scale * 2, 100)
    else:
        xpos = linspace(lognorm.ppf(0.01, shape, scale=scale),
                        lognorm.ppf(0.99, shape, scale=scale), 100)
    plt.plot(xpos, lognorm.pdf(xpos, shape, scale=scale))
    plt.hist(rads, normed=True)
    plt.savefig('packing_histogram.png')
    plt.savefig('packing_histogram.pdf')
    if show_plot:
        plt.show()


def read_results():
    """Reads results of packing algorithm."""
    with open("packing.nfo", "r") as fin:
        fin.readline()
        fin.readline()
        por_theory = float(fin.readline().split()[2])
        por_final = float(fin.readline().split()[2])
        print('Theoretical porosity:', por_theory)
        print('Final porosity:', por_final)
    with open("packing.xyzd", "rb") as fin:
        btxt = fin.read()
        txt = list(struct.unpack("<" + "d" * (len(btxt) // 8), btxt))
        data = array(zip(*[iter(txt)] * 4))
    data[:, 3] = data[:, 3] * \
        ((1 - por_final) / (1 - por_theory))**(1 / 3)
    return data


def render_packing(data, domain=1.0, pixels=1000):
    """Save picture of packed domain. Uses spack.
    https://pyspack.readthedocs.io/en/latest/"""
    pack = spack.Packing(data[:, 0:3], data[:, 3], L=domain)
    print(pack.contacts())
    scene = pack.scene(rot=pi / 4, camera_height=0.5,
                       camera_dist=2.5e1, angle=4, cmap='autumn',
                       floater_color=None)
    scene.render('packing.png', width=pixels,
                 height=pixels, antialiasing=0.0001)


def generate_structure(flag):
    """Runs the packing algorithm."""
    if os.path.isfile("packing.nfo"):
        os.remove(os.path.abspath("packing.nfo"))
    proc = subprocess.Popen(['PackingGeneration.exe', flag])
    proc.wait()
    if not os.path.isfile("packing.nfo"):
        print('Try to change number of particles or size distribution.')
        raise Exception('Packing algorithm failed.')


def pack_spheres(shape, scale, number_of_cells, algorithm):
    """Packs spheres into a periodic domain. Creates Project01.rco with sphere
    centers and radii. Simple model is implemented directly, other algorithms
    use Vasili Baranov's code:
    https://github.com/VasiliBaranov/packing-generation."""
    if algorithm == 'simple':
        data = simple_packing(shape, scale, number_of_cells)
    else:
        create_input(number_of_cells)
        make_csd(shape, scale, number_of_cells)
        generate_structure(algorithm)
        data = read_results()
    np.savetxt('Project01.rco', data)
    render_packing(data)
