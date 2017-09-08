#!/usr/bin/env python
"""Python script, which organizes creation of the foam.

@file       run.py
@namespace  FoamConstruction.run
@ingroup    mod_foamConstruction
@author     Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
@details

First, the geometric tessellation is performed so that the resulting foam has
the correct bubble size distribution. Then several mesh conversions are made to
obtain the foam image in desired format. Finally, foam is voxelized to desired
foam density and struts are optionally added.

Usage:
    simulation.py [-h | --help] [-i input_file] [--verbose]

Options:
    -h --help       Show this screen.
    -i input_file   Json file with inputs. Uses default file otherwise.
    --verbose       Print more information.
"""
from __future__ import division, print_function
import os
import json
import datetime
from blessings import Terminal
from docopt import docopt
from scipy.optimize import minimize_scalar as minimize_scalar
import FoamGeometryConstruction_Periodic
import periodicBox
import vtkconv
import geo_tools
# Creates terminal for colour output
TERM = Terminal()


def porOpt(vx):
    """Objective function.

    For finding size of box, which would give desired porosity.
    @param[in] vx Box size
    """
    filename = INPUTS["filename"]
    porosity = INPUTS["structured_grid_options"]["porosity"]
    vx = int(vx)
    vy = vx
    vz = vx
    if os.path.isfile(filename + 'Box.vtk'):
        os.remove(filename + 'Box.vtk')
    if not os.path.isfile(filename + 'Box.ply'):
        raise SystemError(".ply file is missing. Nothing to binarize.")
    os.system(
        "binvox -e -d {0:d} -rotz -rotx -rotz -rotz -t vtk ".format(vx)
        + filename + "Box.ply >binvox.out"
    )
    with open('binvox.out') as data_file:
        for line in data_file:
            if "counted" in line:
                solidVoxel, totalVoxel =\
                    [int(s) for s in line.split() if s.isdigit()]
                eps = 1 - solidVoxel / totalVoxel
                print("dimension: {0:4d}, porosity: {1:f}".format(vx, eps))
                return (eps - porosity)**2


def porfsOpt(x):
    """Objective function.

    For finding size of box, which would give desired porosity and
    strut content.
    @param[in] x Box size
    """
    global DEDGE
    filename = INPUTS["filename"]
    porosity = INPUTS["structured_grid_options"]["porosity"]
    strutContent = INPUTS["structured_grid_options"]["strutContent"]
    vx = int(x)
    vy = vx
    vz = vx
    if os.path.isfile(filename + 'Box.vtk'):
        os.remove(filename + 'Box.vtk')
    if not os.path.isfile(filename + 'Box.ply'):
        raise SystemError(".ply file is missing. Nothing to binarize.")
    os.system(
        "binvox -e -d {0:d} -rotz -rotx -rotz -rotz -t vtk ".format(vx)
        + filename + "Box.ply >binvox.out"
    )
    filenameIn = filename + "Box.vtk"
    filenameOut = filename + "Box-ascii.vtk"
    dx = dy = dz = INPUTS["packing_options"]["domain_size"]
    origin = [dx, dy, dz]
    spacing = [dx / vx, dy / vy, dz / vz]
    vtkconv.main(filenameIn, filenameOut, origin, spacing)
    f = open("foamreconstr.in", "w")
    f.write("0\n")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("{0:f}\n".format(DEDGE))
    f.write("{0:f}\n".format(1 - strutContent * (1 - porosity)))
    f.write("0\n")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("1\n")
    f.write("0\n")
    f.write("1\n")
    f.write("0\n")
    f.write(filename + "BoxStruts\n")
    f.write(filename + "Box-ascii.vtk\n")
    f.write(filename + ".gnu\n")
    f.write("name\n")
    f.write("descriptors.txt" + "\n")
    f.write("parameters.txt" + "\n")
    f.close()
    os.system("./foamreconstr/foamreconstr")
    f = open("descriptors.txt", "r")
    eps = float(f.readline())
    fs = float(f.readline())
    f.close()
    f = open("parameters.txt", "r")
    DEDGE = float(f.readline())
    f.close()
    resid = ((eps - porosity) / porosity)**2
    print("dimension: {0:4d}, porosity: {1:f}".format(
        vx, eps) + ", strut content: {0:f}".format(fs))
    return resid


def periodic_box(filename, render_box):
    """Uses gmsh, vtk, and meshconv to move closed foam to periodic box."""
    dx = dy = dz = INPUTS["packing_options"]["domain_size"]
    # Convert .geo to .stl
    print(
        TERM.yellow +
        "Convert .geo to .stl" +
        TERM.normal
    )
    os.system("gmsh -n -2 -format stl " + filename + ".geo >gmsh.out")
    # Move to periodic box
    print(
        TERM.yellow +
        "Move to periodic box" +
        TERM.normal
    )
    xmin = dx
    ymin = dy
    zmin = dz
    periodicBox.main(
        filename + ".stl", filename + "Box.stl", xmin, ymin, zmin, dx, dy, dz,
        render_box
    )
    # Convert .stl to .ply
    print(
        TERM.yellow +
        "Convert .stl to .ply" +
        TERM.normal
    )
    os.system("meshconv " + filename + "Box.stl -c ply")


def binarize_box():
    """Creates foam with desired porosity and strut content on structured grid."""
    # Binarize and save as .vtk
    if strutContent == 0:
        # Find the size of box, which would give desired porosity
        # This method is not optimal, since the solver doesn't know that the
        # function takes only integer arguments
        print(
            TERM.yellow +
            "Optimizing porosity" +
            TERM.normal
        )
        res = minimize_scalar(
            porOpt, bracket=[100, 120], method='Brent', tol=1e-2
        )
        vx = res.x
        vx = int(vx)
        print('box size: {0:d}'.format(vx))
        vy = vx
        vz = vx
        print(
            TERM.yellow +
            "Creating and saving optimal foam" +
            TERM.normal
        )
        porOpt(vx)  # Call it with the optimized box size
        # Convert binary .vtk to ascii .vtk
        print(
            TERM.yellow +
            "Convert binary .vtk to ascii .vtk" +
            TERM.normal
        )
        origin = [dx, dy, dz]
        spacing = [dx / vx, dy / vy, dz / vz]
        vtkconv.main(filename + "Box.vtk", filename +
                     "Box-ascii.vtk", origin, spacing)
    else:
        print(
            TERM.yellow +
            "Optimizing porosity and strut content" +
            TERM.normal
        )
        res = minimize_scalar(
            porfsOpt, bracket=[150, 200], method='Brent', tol=1e-2
        )
        # res=minimize_scalar(
        #     porfsOpt,bounds=[200,250],method='bounded',tol=2e0
        # )
        vx = res.x
        vx = int(vx)
        print('optimal box size: {0:d}'.format(vx))
        vy = vx
        vz = vx
        print(
            TERM.yellow +
            "Creating and saving optimal foam" +
            TERM.normal
        )
        if os.path.isfile(filename + 'Box.vtk'):
            os.remove(filename + 'Box.vtk')
        if not os.path.isfile(filename + 'Box.ply'):
            raise SystemError(".ply file is missing. Nothing to binarize.")
        os.system(
            "binvox -e -d {0:d}".format(vx) + " -rotz -rotx -rotz -rotz "
            + "-t vtk " + filename + "Box.ply >binvox.out"
        )
        origin = [dx, dy, dz]
        spacing = [dx / vx, dy / vy, dz / vz]
        vtkconv.main(filename + "Box.vtk", filename +
                     "Box-ascii.vtk", origin, spacing)
        f = open("foamreconstr.in", "w")
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write("{0:f}\n".format(DEDGE))
        f.write("{0:f}\n".format(1 - strutContent * (1 - porosity)))
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write("0\n")
        f.write("0\n")
        f.write("0\n")
        f.write("0\n")
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write(filename + "BoxStruts\n")
        f.write(filename + "Box-ascii.vtk\n")
        f.write(filename + ".gnu\n")
        f.write("name\n")
        f.write("descriptors.txt" + "\n")
        f.write("parameters.txt" + "\n")
        f.close()
        os.system("./foamreconstr/foamreconstr")


def structured_grid():
    """Creates foam discretized on structured grid."""
    if INPUTS["extract_center_cells"]:
        print(
            TERM.yellow +
            "Extracting center cells from tessellation." +
            TERM.normal
        )
        geo_tools.extract_center_cells(
            INPUTS["filename"],
            INPUTS["packing_options"]["number_of_cells"])
    if INPUTS["move_to_periodic_box"]:
        print(
            TERM.yellow +
            "Extracting center cells from tessellation." +
            TERM.normal
        )
        periodic_box(
            INPUTS["filename"],
            INPUTS["structured_grid_options"]["render_box"])


def main():
    """Main function.

    Executed when running the script from command line.
    """
    time_start = datetime.datetime.now()
    dx = dy = dz = INPUTS["packing_options"]["domain_size"]
    if INPUTS["packing"]:
        print(
            TERM.yellow +
            "Packing spheres." +
            TERM.normal
        )
        FoamGeometryConstruction_Periodic.pack_spheres(
            INPUTS["packing_options"]["average_radius"],
            INPUTS["packing_options"]["deviance"],
            INPUTS["packing_options"]["number_of_cells"],
            INPUTS["packing_options"]["alternative_algorithm"])
    if INPUTS["tessellation"]:
        print(
            TERM.yellow +
            "Tessellating." +
            TERM.normal
        )
        FoamGeometryConstruction_Periodic.tessellate(
            INPUTS["filename"],
            INPUTS["packing_options"]["number_of_cells"],
            INPUTS["tessellation_options"]["visualize_tesselation"])
    if INPUTS["structured_grid"]:
        print(
            TERM.yellow +
            "Creating structured grid." +
            TERM.normal
        )
        structured_grid()
    time_end = datetime.datetime.now()
    print("Foam created in: {}".format(time_end - time_start))


if __name__ == "__main__":
    ARGS = docopt(__doc__)
    if ARGS['-i']:
        INPUT_FILE = ARGS['-i']
    else:
        INPUT_FILE = 'input.json'
    with open(INPUT_FILE, 'r') as ifl:
        INPUTS = json.load(ifl)
        DEDGE = INPUTS["structured_grid_options"]["strut_size_guess"]
    main()
