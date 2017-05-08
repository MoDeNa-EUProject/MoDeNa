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

def periodic_surfaces(points, lines, line_loops, surfaces, vec):
    """Returns list of periodic surface pairs."""
    surface_points = [[]]*len(line_loops)
    boundary_points_ind = []
    for surface in surfaces:
        for line in line_loops[surface - 1]:
            for point in lines[line - 1]:
                if point not in surface_points[surface - 1]:
                    surface_points[surface - 1] = surface_points[surface - 1] \
                        + [point]
                if point not in boundary_points_ind:
                    boundary_points_ind.append(point)
    boundary_points = []
    for point in boundary_points_ind:
        boundary_points.append(points[point - 1])
    for i, point in enumerate(surface_points):
        point.sort()
        surface_points[i] = point
    # print(surface_points)
    # print(boundary_points_ind)
    # print(boundary_points)
    eps = 1e-8
    periodic_points = [None]*len(points)
    for i, firstpoint in enumerate(boundary_points):
        for j, secondpoint in enumerate(boundary_points):
            if np.sum(np.abs(firstpoint + vec - secondpoint)) < eps:
                periodic_points[boundary_points_ind[i] - 1] = \
                    boundary_points_ind[j]
    # print(periodic_points)
    psurfs = []
    for i, surf in enumerate(surface_points):
        if surf:
            per_surf = []
            for point in surf:
                per_surf.append(periodic_points[point - 1])
            if None not in per_surf:
                per_surf.sort()
                if per_surf in surface_points:
                    psurfs.append(
                        [i + 1, surface_points.index(per_surf) + 1]
                    )
    # print(psurfs)
    return psurfs

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
    psurfX = periodic_surfaces(
        points, lines, line_loops, surf, np.array([1, 0, 0])
    )
    print(psurfX)
    psurfY = periodic_surfaces(
        points, lines, line_loops, surf, np.array([0, 1, 0])
    )
    print(psurfY)

if __name__ == "__main__":
    main()
