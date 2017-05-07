#!/usr/bin/env python
"""
@brief      Manipulates .geo input files for gmsh.
@author     Pavel Ferkl
"""
from __future__ import print_function
import re
import numpy as np
def my_find_all(regex, text):
    """My definition of findall. Returns top level group in list."""
    matches = re.finditer(regex, text)
    my_list = []
    for match in matches:
        my_list.append(match.group(0))
    return my_list

def read_geo(geo_file):
    """Reads geometry input file for gmsh."""
    with open(geo_file, "r") as text_file:
        text = text_file.read()
        point_s = my_find_all(
            r'Point\s[(][0-9]+[)]\s[=]\s[{]'
            + r'[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)[,]'
            + r'[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)[,]'
            + r'[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)[}][;]',
            text
        )
        line_s = my_find_all(
            r'Line\s[(][0-9]+[)]\s[=]\s[{][0-9]+[,][0-9]+[}][;]', text
        )
        line_loop_s = my_find_all(
            r'Line\sLoop\s[(][0-9]+[)]\s[=]\s[{]([+-]?[0-9]+[,]?)+[}][;]', text
        )
        surface_s = my_find_all(
            r'(Surface\s[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?)+[}][;])'
            + r'(?!.*Physical.*)',
            text
        )
        physical_surface_s = my_find_all(
            r'Physical\sSurface\s[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?)+[}][;]',
            text
        )
        surface_loop_s = my_find_all(
            r'Surface\sLoop\s[(][0-9]+[)]\s[=]\s[{]([+-]?[0-9]+[,]?)+[}][;]',
            text
        )
        volume_s = my_find_all(
            r'Volume\s[(][0-9]+[)]\s[=]\s[{]([0-9]+[,]?)+[}][;]',
            text
        )
        return point_s, line_s, line_loop_s, surface_s,\
            physical_surface_s, surface_loop_s, volume_s

def fix_strings(strings):
    """
    Removes negative sign (orientation); opencascade has problems otherwise.
    """
    for i, line in enumerate(strings):
        strings[i] = re.sub('[-]', '', line)

def save_geo(
        geo_file, point_s, line_s, line_loop_s, surface_s, physical_surface_s,
        surface_loop_s, volume_s, opencascade=True, prepend_plane=False
    ):
    """Saves geometry input file for gmsh."""
    with open(geo_file, "w") as text_file:
        if prepend_plane:
            prepend = "Plane "
        else:
            prepend = ""
        if opencascade:
            text_file.write('SetFactory("OpenCASCADE");\n')
        for line in point_s:
            text_file.write("{}\n".format(line))
        for line in line_s:
            text_file.write("{}\n".format(line))
        for line in line_loop_s:
            text_file.write("{}\n".format(line))
        for line in surface_s:
            text_file.write(prepend + "{}\n".format(line))
        for line in physical_surface_s:
            text_file.write("{}\n".format(line))
        for line in surface_loop_s:
            text_file.write("{}\n".format(line))
        for line in volume_s:
            text_file.write("{}\n".format(line))

def extract(text_array, argtype):
    """
    Returns list of points, lines, or other object. Works for int or float
    variables.
    """
    index = [None]*len(text_array)
    for line in text_array:
        part = line.split("(")
        i = int(part[1].split(")")[0])
        fraction = line.split("{")
        fraction = fraction[1].split("}")
        fraction = fraction[0].split(",")
        fraction = np.array(fraction)
        if argtype == "int":
            fraction = np.absolute(fraction.astype(np.int)).tolist()
        elif argtype == "float":
            fraction = fraction.astype(np.float)
        index[i-1] = fraction
    return index

def surfaces_in_z_plane(points, lines, line_loops, z_coord):
    """Finds surafces that lie completely lie in z plane"""
    points_in_plane = []
    for i, point in enumerate(points):
        if point[2] == z_coord:
            points_in_plane.append(i+1)
    lines_in_plane = []
    for i, line in enumerate(lines):
        if line[0] in points_in_plane and line[1] in points_in_plane:
            lines_in_plane.append(i+1)
    line_loops_in_plane = []
    for i, line_loop in enumerate(line_loops):
        log = True
        for line in line_loop:
            if line not in lines_in_plane:
                log = False
        if log:
            line_loops_in_plane.append(i+1)
    surfaces_in_plane = line_loops_in_plane
    # print(points_in_plane)
    # print(points[points_in_plane[0]])
    # print(lines_in_plane)
    # print(lines[lines_in_plane[0]])
    # print(line_loops_in_plane)
    return surfaces_in_plane

def other_surfaces(surface_loops, surf0, surf1):
    """
    Returns list of boundary surfaces, which are not in surf0 or surf1.
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

def periodic_surafces(points, lines, line_loops, surfaces, vec):
    """Returns list of periodic surface doubles."""
    surface_points = []
    for surface in surfaces:
        for line in line_loops[surface - 1]:
            for point in lines[line - 1]:
                print(points[point - 1])
                exit()

def main():
    """Main subroutine. Just for testing of functionality."""
    point_s, line_s, line_loop_s, surface_s,\
        physical_surface_s, surface_loop_s, volume_s \
        = read_geo("Foam.geo")
    fix_strings(line_loop_s)
    fix_strings(surface_loop_s)
    # save_geo("RVE27_fix.geo", point_s, line_s, line_loop_s, surface_s,\
    #     physical_surface_s, surface_loop_s, volume_s)
    points = extract(point_s, 'float')
    lines = extract(line_s, 'int')
    line_loops = extract(line_loop_s, 'int')
    surfaces = extract(surface_s, 'int')
    physical_surfaces = extract(physical_surface_s, 'int')
    surface_loops = extract(surface_loop_s, 'int')
    volumes = extract(volume_s, 'int')
    surf0 = surfaces_in_z_plane(points, lines, line_loops, 0.0)
    print(surf0)
    surf1 = surfaces_in_z_plane(points, lines, line_loops, 1.0)
    print(surf1)
    surf = other_surfaces(surface_loops, surf0, surf1)
    print(surf)
    periodic_surafces(points, lines, line_loops, surf, [1, 0, 0])

if __name__ == "__main__":
    main()
