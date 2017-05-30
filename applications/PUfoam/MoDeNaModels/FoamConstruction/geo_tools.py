#!/usr/bin/env python
"""
@brief      Manipulates .geo input files for gmsh.
@author     Pavel Ferkl
"""
from __future__ import print_function
import re
import shutil
import subprocess as sp
import numpy as np
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
    """Saves geometry input file for gmsh."""
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
    """Extracts geo data to lists from list of geo strings."""
    edat = {}
    for key in sdat:
        lines = dict()
        for line in sdat[key]:
            part = line.split("(")
            ind = int(part[1].split(")")[0])
            fraction = line.split("{")
            fraction = fraction[1].split("}")
            fraction = fraction[0].split(",")
            if key == "point":
                fraction = np.array(fraction[0:3])
                fraction = fraction.astype(np.float)
                for j, number in enumerate(fraction):
                    if abs(number) < 1e-8:
                        fraction[j] = 0
            else:
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
            for i, j in enumerate(edat[key]):
                sdat[key].append(
                    '{0} {{{1}}} = {{{2}}} Translate{{-1,0,0}};'.format(
                        NAMES[key], j[0], j[1]
                    )
                )
        elif key == 'periodic_surface_Y':
            for i, j in enumerate(edat[key]):
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

def other_surfaces(surface_loops, surf0, surf1):
    """
    Returns list of boundary surfaces, which are not in surf0 or surf1.
    Assumes that inner surfaces are shared by two volumes.
    """
    all_surfaces = []
    for surface_loop in surface_loops:
        all_surfaces += surface_loop
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

def periodic_surfaces(edat, surfaces, vec):
    """Returns list of periodic surface pairs."""
    surface_points = dict()
    boundary_points = dict()
    for surface in surfaces:
        surface_points[surface] = []
        for line in edat['line_loop'][surface]:
            for point in edat['line'][line]:
                if point not in surface_points[surface]:
                    surface_points[surface] = surface_points[surface] + [point]
                if point not in boundary_points:
                    boundary_points[point] = edat['point'][point]
    # boundary_points = []
    # for point in boundary_points_ind:
    #     boundary_points.append(edat['point'][point])
    for i, point in surface_points.iteritems():
        point.sort()
        surface_points[i] = point
    # print(surface_points)
    # print(boundary_points_ind)
    # print(boundary_points)
    # exit()
    eps = 1e-8
    periodic_points = dict()
    for i, firstpoint in boundary_points.iteritems():
        for j, secondpoint in boundary_points.iteritems():
            if np.sum(np.abs(firstpoint + vec - secondpoint)) < eps:
                periodic_points[i] = j
    # print(periodic_points)
    # exit()
    psurfs = []
    for i, surf in surface_points.iteritems():
        per_surf = []
        for point in surf:
            if point in periodic_points:
                per_surf.append(periodic_points[point])
            else:
                per_surf.append(None)
        if None not in per_surf:
            per_surf.sort()
            # print(surf, per_surf)
            if per_surf in surface_points.values():
                psurfs.append(
                    [i, list(surface_points.keys())[list(surface_points.values()).index(per_surf)]]
                )
    # print(psurfs)
    return psurfs

def remove_duplicity(edat):
    """Removes duplicit points, lines, etc."""
    eps = 1e-8
    dupl = dict()
    for i, item1 in edat['point'].iteritems():
        for j, item2 in edat['point'].iteritems():
            if i != j and i > j and np.sum(np.abs(item1 - item2)) < eps:
                if i not in dupl:
                    dupl[i] = []
                dupl[i].append(j)
    for i in dupl:
        del edat['point'][i]
    for values in edat['line'].itervalues():
        for j, value in enumerate(values):
            if value in dupl:
                values[j] = min(dupl[value])
    dupl = dict()
    for i, item1 in edat['line'].iteritems():
        for j, item2 in edat['line'].iteritems():
            if i != j and i > j and sorted(item1) == sorted(item2):
                if i not in dupl:
                    dupl[i] = []
                dupl[i].append(j)
    for i in dupl:
        del edat['line'][i]
    for values in edat['line_loop'].itervalues():
        for j, value in enumerate(values):
            if value in dupl:
                values[j] = min(dupl[value])
    dupl = dict()
    for i, item1 in edat['line_loop'].iteritems():
        for j, item2 in edat['line_loop'].iteritems():
            if i != j and i > j and sorted(item1) == sorted(item2):
                if i not in dupl:
                    dupl[i] = []
                dupl[i].append(j)
    for i in dupl:
        del edat['line_loop'][i]
        del edat['surface'][i]
    for values in edat['surface_loop'].itervalues():
        for j, value in enumerate(values):
            if value in dupl:
                values[j] = min(dupl[value])
    # dupl = dict() # no duplicit volumes
    # for i, item1 in edat['surface_loop'].iteritems():
    #     for j, item2 in edat['surface_loop'].iteritems():
    #         if i != j and i > j and sorted(item1) == sorted(item2):
    #             if i not in dupl:
    #                 dupl[i] = []
    #             dupl[i].append(j)
    # print(dupl)
    # exit()

def move_to_box(infile, wfile, outfile, volumes):
    """Moves periodic closed foam to periodic box."""
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
            + '{{Volume{{1:{0}}};}}'.format(mvol)
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 1)
        )
        wfl.write(
            'zoh() = BooleanIntersection'
            + '{{Volume{{1:{0}}};}}'.format(mvol)
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 2)
        )
        wfl.write(
            'zin() = BooleanIntersection'
            + '{{Volume{{1:{0}}}; Delete;}}'.format(mvol)
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 3)
        )
        wfl.write('Translate{0,0, 1}{Volume{zol()};}\n')
        wfl.write('Translate{0,0,-1}{Volume{zoh()};}\n\n')
        wfl.write(
            'yol() = BooleanIntersection'
            + '{Volume{zol(),zoh(),zin()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 4)
        )
        wfl.write(
            'yoh() = BooleanIntersection'
            + '{Volume{zol(),zoh(),zin()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 5)
        )
        wfl.write(
            'yin() = BooleanIntersection'
            + '{Volume{zol(),zoh(),zin()}; Delete;}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 6)
        )
        wfl.write('Translate{0, 1,0}{Volume{yol()};}\n')
        wfl.write('Translate{0,-1,0}{Volume{yoh()};}\n\n')
        wfl.write(
            'xol() = BooleanIntersection'
            + '{Volume{yol(),yoh(),yin()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 7)
        )
        wfl.write(
            'xoh() = BooleanIntersection'
            + '{Volume{yol(),yoh(),yin()};}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 8)
        )
        wfl.write(
            'xin() = BooleanIntersection'
            + '{Volume{yol(),yoh(),yin()}; Delete;}'
            + '{{Volume{{{0}}}; Delete;}};\n'.format(mvol + 9)
        )
        wfl.write('Translate{ 1,0,0}{Volume{xol()};}\n')
        wfl.write('Translate{-1,0,0}{Volume{xoh()};}\n\n')
    call = sp.Popen(['gmsh', wfile, '-0'])
    call.wait()
    shutil.move(wfile+'_unrolled', outfile)

def main():
    """Main subroutine. Just for testing of functionality."""
    # sdat = read_geo("FoamClosed.geo") # string data
    # fix_strings(sdat['line_loop'])
    # fix_strings(sdat['surface_loop'])
    # sdat.pop('physical_surface')
    # save_geo("FoamClosedFixed.geo", sdat)
    # move_to_box(
    #     "FoamClosedFixed.geo", "move_to_box.geo", "FoamBox.geo",
    #     range(1, len(sdat['volume']) + 1)
    # )
    sdat = read_geo("FoamBox.geo") # string data
    edat = extract_data(sdat) # extracted data
    remove_duplicity(edat)
    surf0 = surfaces_in_plane(edat, 0.0, 2)
    print(surf0)
    surf1 = surfaces_in_plane(edat, 1.0, 2)
    print(surf1)
    # surf = other_surfaces(edat['surface_loop'].itervalues(), surf0, surf1)
    surf = surfaces_in_plane(edat, 0.0, 1) + surfaces_in_plane(edat, 1.0, 1) \
        + surfaces_in_plane(edat, 0.0, 0) + surfaces_in_plane(edat, 1.0, 0)
    print(surf)
    edat.pop('physical_surface')
    edat['physical_surface'] = {1:surf0, 2:surf1, 3:surf}
    edat['periodic_surface_X'] = periodic_surfaces(
        edat, surf, np.array([1, 0, 0])
    )
    print(edat['periodic_surface_X'])
    edat['periodic_surface_Y'] = periodic_surfaces(
        edat, surf, np.array([0, 1, 0])
    )
    print(edat['periodic_surface_Y'])
    edat['physical_volume'] = {1:edat['volume'].keys()}
    # edat['physical_volume'] = edat['volume']
    print(edat.keys())
    sdat2 = collect_strings(edat)
    print(sdat2.keys())
    save_geo("FoamBoxFixed.geo", sdat2)

if __name__ == "__main__":
    main()
