#!/usr/bin/env python
"""
@brief      Manipulates .geo input files for gmsh.
@author     Pavel Ferkl
"""
from __future__ import print_function
import re
def my_find_all(regex, text):
    """My definition of findall. Returns top level group in list."""
    matches = re.finditer(regex, text)
    my_list = []
    for match in matches:
        my_list.append(match.group(0))
    return my_list

def read_geo(geo_file):
    """reads geometry input file for gmsh"""
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
    removes negative sign (orientation); opencascade has problems otherwise
    """
    for i, line in enumerate(strings):
        strings[i] = re.sub('[-]', '', line)

def save_geo(geo_file, point_s, line_s, line_loop_s, surface_s,\
        physical_surface_s, surface_loop_s, volume_s):
    """saves geometry input file for gmsh"""
    with open(geo_file, "w") as text_file:
        for line in point_s:
            text_file.write("{}\n".format(line))
        for line in line_s:
            text_file.write("{}\n".format(line))
        for line in line_loop_s:
            text_file.write("{}\n".format(line))
        for line in surface_s:
            text_file.write("{}\n".format(line))
        for line in physical_surface_s:
            text_file.write("{}\n".format(line))
        for line in surface_loop_s:
            text_file.write("{}\n".format(line))
        for line in volume_s:
            text_file.write("{}\n".format(line))

def main():
    """main subroutine"""
    point_s, line_s, line_loop_s, surface_s,\
        physical_surface_s, surface_loop_s, volume_s \
        = read_geo("RVE27.geo")
    fix_strings(line_loop_s)
    fix_strings(surface_loop_s)
    save_geo("RVE27_fix.geo", point_s, line_s, line_loop_s, surface_s,\
        physical_surface_s, surface_loop_s, volume_s)

if __name__ == "__main__":
    main()
