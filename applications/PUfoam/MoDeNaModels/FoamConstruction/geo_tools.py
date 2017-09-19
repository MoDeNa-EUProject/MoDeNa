"""Manipulates .geo input files for gmsh.

@author Pavel Ferkl
"""
from __future__ import print_function
import re
import shutil
import subprocess as sp
import numpy as np
from blessings import Terminal
from docopt import docopt
NAMES = {
    'point': 'Point',
    'line': 'Line',
    'line_loop': 'Line Loop',
    'surface': 'Plane Surface',
    'surface_loop': 'Surface Loop',
    'volume': 'Volume',
    'periodic_surface_X': 'Periodic Surface',
    'periodic_surface_Y': 'Periodic Surface',
    'physical_surface': 'Physical Surface',
    'physical_volume': 'Physical Volume'
}
NAME_LIST = [
    'point',
    'line',
    'line_loop',
    'surface',
    'surface_loop',
    'volume',
    'periodic_surface_X',
    'periodic_surface_Y',
    'physical_surface',
    'physical_volume'
]


def my_find_all(regex, text):
    """My definition of findall. Returns top level group in list."""
    matches = re.finditer(regex, text)
    my_list = []
    for match in matches:
        my_list.append(match.group(0))
    return my_list


def read_geo(geo_file, ignore_point_format=True, plane_surface=True):
    """
    Reads geometry input file for gmsh into dictionary. Based on regular
    expressions. Points can contain optional fourth argument, thus it is better
    to include everything in curly braces. Some geo files use Surface, some
    Plane Surface. You should specify what you want to read.
    """
    with open(geo_file, "r") as text_file:
        text = text_file.read()
        sdat = {}
        rexp = {}
        if ignore_point_format:
            rexp['point'] = r'Point\s?[(][0-9]+[)]\s[=]\s[{](.*?)[}][;]'
        else:
            rexp['point'] = (
                r'Point\s[(][0-9]+[)]\s[=]\s[{]'
                + r'[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)[,]'
                + r'[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)[,]'
                + r'[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)[}][;]'
            )
        rexp['line'] = r'Line\s?[(][0-9]+[)]\s[=]\s[{][0-9]+[,]\s?[0-9]+[}][;]'
        rexp['line_loop'] = (
            r'Line\sLoop\s?[(][0-9]+[)]\s[=]\s[{]([+-]?[0-9]+[,]?\s?)+[}][;]'
        )
        if plane_surface:
            rexp['surface'] = (
                r'Plane\sSurface\s?[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?\s?)+[}][;]'
            )
        else:
            rexp['surface'] = (
                r'(Surface\s[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?)+[}][;])'
                + r'(?!.*Physical.*)',
            )
        rexp['physical_surface'] = (
            r'Physical\sSurface\s?[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?\s?)+[}][;]'
        )
        rexp['surface_loop'] = (
            r'Surface\sLoop\s?[(][0-9]+[)]\s[=]\s[{]([+-]?[0-9]+[,]?\s?)+[}][;]'
        )
        rexp['volume'] = (
            r'Volume\s?[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?\s?)+[}][;]'
        )
        rexp['physical_volume'] = (
            r'Physical\sVolume\s?[(]["][a-z]+["][)]\s[=]\s'
            + r'[{]([0-9]+[,]?\s?)+[}][;]'
        )
        for key in rexp:
            sdat[key] = my_find_all(rexp[key], text)
        return sdat


def fix_strings(strings):
    """
    Removes negative sign (orientation) from loops. OpenCASCADE has problems
    otherwise.
    """
    for i, line in enumerate(strings):
        strings[i] = re.sub('[-]', '', line)


def save_geo(geo_file, sdat, opencascade=True, char_length=0.1):
    """
    Creates geometry input file for gmsh. Input is a dictionary with prepared
    string lines.
    """
    with open(geo_file, "w") as text_file:
        if opencascade:
            text_file.write('SetFactory("OpenCASCADE");\n')
            text_file.write(
                'Mesh.CharacteristicLengthMax = {0};\n'.format(char_length)
            )
        for key in NAME_LIST:
            if key in sdat:
                for line in sdat[key]:
                    text_file.write("{}\n".format(line))


def extract_data(sdat):
    """Extracts geo data to dictionaries from list of geo strings."""
    edat = {}
    for key in sdat:
        lines = dict()
        for line in sdat[key]:
            part = line.split("(")
            if key == "physical_volume":
                ind = part[1].split(")")[0]  # ID of the element
                if ind == '"cells"':
                    ind = 1
                elif ind == '"walls"':
                    ind = 2
            else:
                ind = int(part[1].split(")")[0])  # ID of the element
            fraction = line.split("{")
            fraction = fraction[1].split("}")
            fraction = fraction[0].split(",")
            if key == "point":  # point data consists of floats
                # ignore the optional fourth argument (defines mesh coarseness)
                fraction = np.array(fraction[0:3])
                fraction = fraction.astype(np.float)
                for j, number in enumerate(fraction):
                    if abs(number) < 1e-8:
                        fraction[j] = 0
            else:  # other data consists of integers
                fraction = np.array(fraction)
                fraction = np.absolute(fraction.astype(np.int)).tolist()
            lines[ind] = fraction
        edat[key] = lines
    return edat


def collect_strings(edat):
    """Creates lists of geo strings from geo data."""
    sdat = {}
    for key in edat:
        sdat[key] = []
        if key == 'periodic_surface_X':
            for j in edat[key]:
                sdat[key].append(
                    '{0} {{{1}}} = {{{2}}} Translate{{-1,0,0}};'.format(
                        NAMES[key], j[0], j[1]
                    )
                )
        elif key == 'periodic_surface_Y':
            for j in edat[key]:
                sdat[key].append(
                    '{0} {{{1}}} = {{{2}}} Translate{{0,-1,0}};'.format(
                        NAMES[key], j[0], j[1]
                    )
                )
        else:
            for i, j in edat[key].iteritems():
                j = ','.join(str(e) for e in j)
                sdat[key].append('{0} ({1}) = {{{2}}};'.format(
                    NAMES[key], i, j
                ))
    return sdat


def surfaces_in_plane(edat, coord, direction):
    """Finds surfaces that lie completely lie in a plane"""
    points_in_plane = []
    for i, point in edat['point'].iteritems():
        if point[direction] == coord:
            points_in_plane.append(i)
    lines_in_plane = []
    for i, line in edat['line'].iteritems():
        if line[0] in points_in_plane and line[1] in points_in_plane:
            lines_in_plane.append(i)
    line_loops_in_plane = []
    for i, line_loop in edat['line_loop'].iteritems():
        log = True
        for line in line_loop:
            if line not in lines_in_plane:
                log = False
        if log:
            line_loops_in_plane.append(i)
    return line_loops_in_plane


def other_surfaces(edat, surf0, surf1):
    """
    Returns list of boundary surfaces, which are not in surf0 or surf1. Assumes
    that inner surfaces are shared by two volumes. Remove duplicates before
    calling this function.
    """
    all_surfaces = []
    for surface_loops in edat['volume'].itervalues():
        for surface_loop in surface_loops:
            for surfaces in edat['surface_loop'][surface_loop]:
                all_surfaces += edat['surface'][surfaces]
    count = dict()
    for surface in all_surfaces:
        if surface in count:
            count[surface] += 1
        else:
            count[surface] = 1
    surf = [
        i for i, j in count.items() if j == 1
        and i not in surf0 and i not in surf1
    ]
    return surf


def periodic_surfaces(edat, surfaces, vec, eps=1e-8):
    """
    Returns list of periodic surface pairs in specified direction. Parameter eps
    is used to compare coordinates (floating point numbers).
    """
    surface_points = dict()  # point IDs for each boundary surface
    boundary_points = dict()  # dictionary with only boundary points
    for surface in surfaces:
        surface_points[surface] = []
        for line in edat['line_loop'][surface]:
            for point in edat['line'][line]:
                if point not in surface_points[surface]:
                    surface_points[surface] += [point]
                if point not in boundary_points:
                    boundary_points[point] = edat['point'][point]
    # sort point IDs so that you can compare later
    for point in surface_points.itervalues():
        point.sort()
    # dictionary with ID of periodic point for each point that has one
    periodic_points = dict()
    for i, point in boundary_points.iteritems():
        for j, secondpoint in boundary_points.iteritems():
            if np.sum(np.abs(point + vec - secondpoint)) < eps:
                periodic_points[i] = j
    psurfs = []  # list of periodic surface pairs (IDs)
    for i, surface in surface_points.iteritems():
        # Try to create surface using IDs of periodic points. Use None if there
        # is no periodic point in specified direction.
        per_surf = [
            periodic_points[point] if point in periodic_points else None
            for point in surface
        ]
        if None not in per_surf:
            per_surf.sort()  # sort so you can find it
            # use ID of current surface and find ID of periodic surface
            psurfs.append(
                [
                    i,
                    surface_points.keys()[
                        surface_points.values().index(per_surf)
                    ]
                ]
            )
    return psurfs


def identify_duplicity(edat, key, number, eps):
    """
    Core algorithm for removing duplicities. User should call remove_duplicity()
    instead. Parameter eps is used to compare coordinates (floating point
    numbers).
    """
    dupl = dict()
    if number == 'float':
        for i, item1 in edat[key].iteritems():
            for j, item2 in edat[key].iteritems():
                if i != j and i > j and np.sum(np.abs(item1 - item2)) < eps:
                    if i not in dupl:
                        dupl[i] = []
                    dupl[i].append(j)
    elif number == 'integer':
        for i, item1 in edat[key].iteritems():
            for j, item2 in edat[key].iteritems():
                if i != j and i > j and sorted(item1) == sorted(item2):
                    if i not in dupl:
                        dupl[i] = []
                    dupl[i].append(j)
    else:
        raise Exception('number argument must be float or integer')
    return dupl


def remove_duplicit_ids_from_keys(edat, dupl, key):
    """Removes duplicit IDs from IDs of entities."""
    for i in dupl:
        del edat[key][i]


def remove_duplicit_ids_from_values(edat, dupl, key):
    """Removes duplicit IDs from values of entities."""
    for values in edat[key].itervalues():
        for j, value in enumerate(values):
            if value in dupl:
                values[j] = min(dupl[value])


def remove_duplicity(edat, eps=1e-10):
    """
    Removes duplicit points, lines, etc. Parameter eps is used to compare
    coordinates (floating point numbers).
    """
    # points
    dupl = identify_duplicity(edat, 'point', 'float', eps)
    remove_duplicit_ids_from_keys(edat, dupl, 'point')
    remove_duplicit_ids_from_values(edat, dupl, 'line')
    # lines
    dupl = identify_duplicity(edat, 'line', 'integer', eps)
    remove_duplicit_ids_from_keys(edat, dupl, 'line')
    remove_duplicit_ids_from_values(edat, dupl, 'line_loop')
    # line loops
    dupl = identify_duplicity(edat, 'line_loop', 'integer', eps)
    remove_duplicit_ids_from_keys(edat, dupl, 'line_loop')
    remove_duplicit_ids_from_keys(edat, dupl, 'surface')
    remove_duplicit_ids_from_values(edat, dupl, 'surface_loop')
    # there are no duplicit volumes


def split_loops(edat, key):
    """
    Makes sure that line and surface loops contain only one loop. Surfaces and
    volumes with holes are instead defined in Surface and Volume entries,
    respectively. Needed because gmsh unrolls geometry in a way, which is
    unusable with OpenCASCADE kernel.
    """
    if key == 'line_loop':
        key2 = 'surface'
    elif key == 'surface_loop':
        key2 = 'volume'
    else:
        raise Exception('can be called only for line_loop or surface_loop')
    for i, item1 in edat[key].iteritems():
        for j, item2 in edat[key].iteritems():
            if i != j and set(item2).issubset((set(item1))):
                for value in item2:
                    item1.remove(value)
                edat[key][i] = item1
                edat[key2][i] = [i, j]
                break


def move_to_box(infile, wfile, outfile, volumes):
    """
    Moves periodic closed foam to periodic box. Uses gmsh, specifically boolean
    operations and transformations from OpenCASCADE. The result is unrolled to
    another geo file so that it can be quickly read and worked with in the
    follow-up work. Operations are performed two times. First for walls (first
    half of volumes) and then for cells.
    """
    with open(wfile, 'w') as wfl:
        mvol = max(volumes)
        wfl.write('SetFactory("OpenCASCADE");\n\n')
        wfl.write('Include "{0}";\n\n'.format(infile))
        wfl.write('Block({0}) = {{-1,-1,-1,3,3,1}};\n'.format(mvol + 1))
        wfl.write('Block({0}) = {{-1,-1, 1,3,3,1}};\n'.format(mvol + 2))
        wfl.write('Block({0}) = {{-1,-1, 0,3,3,1}};\n'.format(mvol + 3))
        wfl.write('Block({0}) = {{-1,-1,-1,3,1,3}};\n'.format(mvol + 4))
        wfl.write('Block({0}) = {{-1, 1,-1,3,1,3}};\n'.format(mvol + 5))
        wfl.write('Block({0}) = {{-1, 0,-1,3,1,3}};\n'.format(mvol + 6))
        wfl.write('Block({0}) = {{-1,-1,-1,1,3,3}};\n'.format(mvol + 7))
        wfl.write('Block({0}) = {{ 1,-1,-1,1,3,3}};\n'.format(mvol + 8))
        wfl.write('Block({0}) = {{ 0,-1,-1,1,3,3}};\n'.format(mvol + 9))
        wfl.write('\n')
        wfl.write(
            'zol() = BooleanIntersection'
            + '{{Volume{{1:{0}}};}}'.format(mvol / 2)
            + '{{Volume{{{0}}};}};\n'.format(mvol + 1)
        )
        wfl.write(
            'zoh() = BooleanIntersection'
            + '{{Volume{{1:{0}}};}}'.format(mvol / 2)
            + '{{Volume{{{0}}};}};\n'.format(mvol + 2)
        )
        wfl.write(
            'zin() = BooleanIntersection'
            + '{{Volume{{1:{0}}}; Delete;}}'.format(mvol / 2)
            + '{{Volume{{{0}}};}};\n'.format(mvol + 3)
        )
        wfl.write('Translate{0,0, 1}{Volume{zol()};}\n')
        wfl.write('Translate{0,0,-1}{Volume{zoh()};}\n\n')
        wfl.write(
            'yol() = BooleanIntersection'
            + '{Volume{zol(),zoh(),zin()};}'
            + '{{Volume{{{0}}};}};\n'.format(mvol + 4)
        )
        wfl.write(
            'yoh() = BooleanIntersection'
            + '{Volume{zol(),zoh(),zin()};}'
            + '{{Volume{{{0}}};}};\n'.format(mvol + 5)
        )
        wfl.write(
            'yin() = BooleanIntersection'
            + '{Volume{zol(),zoh(),zin()}; Delete;}'
            + '{{Volume{{{0}}};}};\n'.format(mvol + 6)
        )
        wfl.write('Translate{0, 1,0}{Volume{yol()};}\n')
        wfl.write('Translate{0,-1,0}{Volume{yoh()};}\n\n')
        wfl.write(
            'xol() = BooleanIntersection'
            + '{Volume{yol(),yoh(),yin()};}'
            + '{{Volume{{{0}}};}};\n'.format(mvol + 7)
        )
        wfl.write(
            'xoh() = BooleanIntersection'
            + '{Volume{yol(),yoh(),yin()};}'
            + '{{Volume{{{0}}};}};\n'.format(mvol + 8)
        )
        wfl.write(
            'xin() = BooleanIntersection'
            + '{Volume{yol(),yoh(),yin()}; Delete;}'
            + '{{Volume{{{0}}};}};\n'.format(mvol + 9)
        )
        wfl.write('Translate{ 1,0,0}{Volume{xol()};}\n')
        wfl.write('Translate{-1,0,0}{Volume{xoh()};}\n\n')
        wfl.write(
            'zol2() = BooleanIntersection'
            + '{{Volume{{{0}:{1}}};}}'.format(mvol / 2 + 1, mvol)
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 1)
        )
        wfl.write(
            'zoh2() = BooleanIntersection'
            + '{{Volume{{{0}:{1}}};}}'.format(mvol / 2 + 1, mvol)
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 2)
        )
        wfl.write(
            'zin2() = BooleanIntersection'
            + '{{Volume{{{0}:{1}}}; Delete;}}'.format(mvol / 2 + 1, mvol)
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 3)
        )
        wfl.write('Translate{0,0, 1}{Volume{zol2()};}\n')
        wfl.write('Translate{0,0,-1}{Volume{zoh2()};}\n\n')
        wfl.write(
            'yol2() = BooleanIntersection'
            + '{Volume{zol2(),zoh2(),zin2()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 4)
        )
        wfl.write(
            'yoh2() = BooleanIntersection'
            + '{Volume{zol2(),zoh2(),zin2()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 5)
        )
        wfl.write(
            'yin2() = BooleanIntersection'
            + '{Volume{zol2(),zoh2(),zin2()}; Delete;}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 6)
        )
        wfl.write('Translate{0, 1,0}{Volume{yol2()};}\n')
        wfl.write('Translate{0,-1,0}{Volume{yoh2()};}\n\n')
        wfl.write(
            'xol2() = BooleanIntersection'
            + '{Volume{yol2(),yoh2(),yin2()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 7)
        )
        wfl.write(
            'xoh2() = BooleanIntersection'
            + '{Volume{yol2(),yoh2(),yin2()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 8)
        )
        wfl.write(
            'xin2() = BooleanIntersection'
            + '{Volume{yol2(),yoh2(),yin2()}; Delete;}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 9)
        )
        wfl.write('Translate{ 1,0,0}{Volume{xol2()};}\n')
        wfl.write('Translate{-1,0,0}{Volume{xoh2()};}\n\n')
        wfl.write('Physical Volume ("walls") = {xol(),xoh(),xin()};\n')
        wfl.write('Physical Volume ("cells") = {xol2(),xoh2(),xin2()};\n\n')
    call = sp.Popen(['gmsh', wfile, '-0'])
    call.wait()
    shutil.move(wfile + '_unrolled', outfile)


def create_walls(edat, wall_thickness=0.01):
    """Creates walls."""
    volume_points = dict()  # point IDs for each volume
    for volume in edat['surface_loop']:
        volume_points[volume] = []
        for surface in edat['surface_loop'][volume]:
            for line in edat['line_loop'][surface]:
                for point in edat['line'][line]:
                    if point not in volume_points[volume]:
                        volume_points[volume] += [point]
        volume_points[volume].sort()
    centroids = dict()  # centroid for each volume
    for volume in edat['surface_loop']:
        total = 0
        for point in volume_points[volume]:
            total += edat['point'][point]
        total /= len(volume_points[volume])
        centroids[volume] = total
    npoints = len(edat['point'])
    nlines = len(edat['line'])
    nsurfaces = len(edat['line_loop'])
    nvolumes = len(edat['surface_loop'])
    for volume in edat['surface_loop'].keys():
        point_map = dict()  # mapping of old points to new points
        nvolumes += 1
        edat['surface_loop'][nvolumes] = []
        for point in volume_points[volume]:
            npoints += 1
            edat['point'][npoints] = edat['point'][point] + wall_thickness * (
                centroids[volume] - edat['point'][point])
            point_map[point] = npoints
        for surface in edat['surface_loop'][volume]:
            nsurfaces += 1
            edat['line_loop'][nsurfaces] = []
            for line in edat['line_loop'][surface]:
                nlines += 1
                edat['line'][nlines] = [
                    point_map[edat['line'][line][0]],
                    point_map[edat['line'][line][1]],
                ]
                edat['line_loop'][nsurfaces] += [nlines]
            edat['surface'][nsurfaces] = [nsurfaces]
            edat['surface_loop'][nvolumes] += [nsurfaces]
        edat['volume'][nvolumes] = [nvolumes]
        edat['volume'][volume] += [nvolumes]
    remove_duplicity(edat)


def extract_center_cells(filename, number_of_cells):
    """Extracts cells in the center from 27 times larger tessellated doamin."""
    sdat = read_geo("{0}RVE27.geo".format(filename))
    edat = extract_data(sdat)
    edges = edat['line'].values()
    faces = edat['line_loop'].values()
    volumes = edat['surface_loop'].values()
    #####################################################
    max0 = list(range(0, number_of_cells))
    for i in range(number_of_cells):
        max0[i] = max(volumes[i])
    max_index_of_faces = max(max0)
    max1 = list(range(0, max_index_of_faces))
    for i in range(max_index_of_faces):
        max1[i] = max(faces[i])
    max_index_of_edges = max(max1)
    max2 = list(range(0, max_index_of_edges))
    for i in range(max_index_of_edges):
        max2[i] = max(edges[i])
    max_index_of_nodes = max(max2)
    ####################################################
    # Making GEO file containing Periodic RVE
    sdat['point'] = sdat['point'][:max_index_of_nodes]
    sdat['line'] = sdat['line'][:max_index_of_edges]
    sdat['line_loop'] = sdat['line_loop'][:max_index_of_faces]
    sdat['surface'] = sdat['surface'][:max_index_of_faces]
    sdat['physical_surface'] = sdat['physical_surface'][:max_index_of_faces]
    sdat['surface_loop'] = sdat['surface_loop'][:number_of_cells]
    sdat['volume'] = sdat['volume'][:number_of_cells]
    save_geo(
        "{0}.geo".format(filename),
        sdat,
        opencascade=False
    )


def main(fname, wall_thickness, verbose):
    """
    Main subroutine. Organizes workflow.

    File.geo -> FileFixed.geo -> FileBox.geo -> FileBoxFixed.geo
    """
    term = Terminal()
    print(
        term.yellow
        + "Working on file {}.geo.".format(fname)
        + term.normal
    )
    # read Neper foam
    sdat = read_geo(fname + ".geo")  # string data
    # Neper creates physical surfaces, which we don't want
    sdat.pop('physical_surface')
    # remove orientation, OpenCASCADE compatibility
    fix_strings(sdat['line_loop'])
    fix_strings(sdat['surface_loop'])
    # create walls
    edat = extract_data(sdat)
    create_walls(edat, wall_thickness)
    sdat = collect_strings(edat)
    save_geo(fname + "Walls.geo", sdat)
    # move foam to a periodic box and save it to a file
    move_to_box(
        fname + "Walls.geo", "move_to_box.geo", fname + "WallsBox.geo",
        range(1, len(sdat['volume']) + 1)
    )
    # read boxed foam
    sdat = read_geo(fname + "WallsBox.geo")  # string data
    edat = extract_data(sdat)  # extracted data
    # duplicity of points, lines, etc. was created during moving to a box
    remove_duplicity(edat)
    # restore OpenCASCADE compatibility
    split_loops(edat, 'line_loop')
    split_loops(edat, 'surface_loop')
    # identification of physical surfaces for boundary conditions
    surf0 = surfaces_in_plane(edat, 0.0, 2)
    if verbose:
        print('Z=0 surface IDs: {}'.format(surf0))
    surf1 = surfaces_in_plane(edat, 1.0, 2)
    if verbose:
        print('Z=1 surface IDs: {}'.format(surf1))
    surf = other_surfaces(edat, surf0, surf1)
    if verbose:
        print('other boundary surface IDs: {}'.format(surf))
    """
    Physical surfaces create problems in mesh conversion step. Bug in gmsh?
    Boundaries will be defined in fenics/dolfin directly.
    TODO: fix this
    """
    # edat['physical_surface'] = {1:surf0, 2:surf1, 3:surf}
    # identification of periodic surfaces for periodic mesh creation
    edat['periodic_surface_X'] = periodic_surfaces(
        edat, surf, np.array([1, 0, 0])
    )
    if verbose:
        print(
            'surface IDs periodic in X: {}'.format(edat['periodic_surface_X'])
        )
    edat['periodic_surface_Y'] = periodic_surfaces(
        edat, surf, np.array([0, 1, 0])
    )
    if verbose:
        print(
            'surface IDs periodic in Y: {}'.format(edat['periodic_surface_Y'])
        )
    # save the final foam
    sdat = collect_strings(edat)
    save_geo(fname + "WallsBoxFixed.geo", sdat)
    print(
        term.yellow
        + "Prepared file {}WallsBoxFixed.geo.".format(fname)
        + term.normal
    )
