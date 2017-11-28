#!/usr/bin/env python
"""Creates input file, so that foamConductivity can be initialized."""
from __future__ import print_function
import sys
import os
import json
from numpy import loadtxt, eye
from modena import SurrogateModel, IndexSet


def my_help():
    "print usage of the module"
    print('Usage: prep_init_foamConductivity v')
    print('v: integer')
    print("0 - try to initialize with points so that foamConductivity " +
          "does not go out of range")
    print('1 - initialize with points from foamAging.json')
    print('2 - initialize with points from foamAging.json except for gas '
          + 'composition - foamConductivity should not go out of range in '
          + 'foamAging simulation')
    return


def no_refitting():
    """
    initialize with dummy data that should ensure that foamConductivity
    does not go out of range
    """
    # this can fail when we change species in gas conductivity importing gas
    # conductivity surrogate model and running this script again should work
    try:
        indxs = list(IndexSet.load("gas_thermal_conductivity_species").names)
    except:
        indxs = ['CO2', 'CyP', 'O2', 'N2', 'Opt', 'Sol']
    ini = {
        "eps": [0.9, 0.0, 0.96, 0.99, 0.7, 0.5],
        "dcell": [200e-6, 0.0, 300e-6, 100e-6, 1e-2, 200e-6],
        "fstrut": [0.0, 1.0, 0.7, 0.6, 0.0, 0.9],
        "T": [280, 549, 300, 350, 330, 300]
    }
    mat = eye(6)
    for i, indx in enumerate(indxs):
        ini["x[" + indx + "]"] = mat[i % 6, :].tolist()
    with open("./inputs/init_foamConductivity.json", "w") as fl:
        json.dump(ini, fl, indent=4, sort_keys=True)
    return


def setIP(a0):
    "creates a list with small perturbation around the input value"
    a = []
    for i in xrange(4):
        a.append(a0)
    upar = 1 - 1e-4
    opar = 1 + 1e-4
    a[0] = a[0] * upar
    if a0 == 0:
        a[1] = 1e-4
    else:
        a[1] = a[1] * opar
    return a


def setIP_fractions(a0):
    "creates a list with small perturbation around the input value"
    a = []
    for i in xrange(4):
        a.append(a0)
    upar = 1 - 1e-4
    opar = 1 + 1e-4
    if a0 > 1:
        a0 = 1.0
    a[0] = a[0] * upar
    if a0 == 0:
        a[1] = 1e-4
    else:
        a[1] = a[1] * opar
        if a[1] > 1:
            a[1] = 1.0
    return a


def get_rho(inputs, results):
    """Extract foam density."""
    if inputs['sourceOfProperty']['foamDensity'] == 'DirectInput':
        rho0 = inputs['morphology']['foamDensity']
    elif inputs['sourceOfProperty']['foamDensity'] == 'BubbleGrowth':
        with open(results + "/bubbleGrowth/after_foaming.txt") as fl2:
            rho0, dcell, xCyP, xCO2 = loadtxt(fl2, skiprows=1, unpack=True)
    elif inputs['sourceOfProperty']['foamDensity'] == 'Qmom0D':
        with open(results + "/CFD0D/after_foaming.txt") as fl2:
            rho0, dcell, xCyP, xCO2 = loadtxt(fl2, skiprows=1, unpack=True)
    elif inputs['sourceOfProperty']['foamDensity'] == 'Qmom3D':
        with open(results + "/CFD3D/after_foaming.txt") as fl2:
            rho0, dcell, xCyP, xCO2 = loadtxt(fl2, skiprows=1, unpack=True)
    else:
        raise Exception("unknown source for foam density")
    return rho0


def get_dcell(inputs, results):
    """Extract cell size."""
    if inputs['sourceOfProperty']['cellSize'] == 'DirectInput':
        dcell0 = inputs['morphology']['cellSize']
    elif inputs['sourceOfProperty']['cellSize'] == 'BubbleGrowth':
        with open(results + "/bubbleGrowth/after_foaming.txt") as fl2:
            rho, dcell0, xCyP, xCO2 = loadtxt(fl2, skiprows=1, unpack=True)
    elif inputs['sourceOfProperty']['cellSize'] == 'Qmom0D':
        with open(results + "/CFD0D/after_foaming.txt") as fl2:
            rho, dcell0, xCyP, xCO2 = loadtxt(fl2, skiprows=1, unpack=True)
    elif inputs['sourceOfProperty']['cellSize'] == 'Qmom3D':
        with open(results + "/CFD3D/after_foaming.txt") as fl2:
            rho, dcell0, xCyP, xCO2 = loadtxt(fl2, skiprows=1, unpack=True)
    else:
        raise Exception("unknown source for cell size")
    return dcell0


def get_gas_comp(inputs, results):
    """Extract gas composition."""
    names = inputs['foamCondition']['initialComposition'].keys()
    if inputs['sourceOfProperty']['gasComposition'] == 'DirectInput':
        xfrac = {}
        for name in names:
            xfrac[name] = inputs['foamCondition']['initialComposition'][name]
    elif inputs['sourceOfProperty']['gasComposition'] == 'BubbleGrowth':
        with open(results + "/bubbleGrowth/after_foaming.txt") as fl2:
            rho, dcell, xCyP0, xCO20 = loadtxt(fl2, skiprows=1, unpack=True)
        for name in names:
            xfrac[name] = 0
        xfrac['CO2'] = xCO20
        xfrac['CyP'] = xCyP0
    elif inputs['sourceOfProperty']['gasComposition'] == 'Qmom0D':
        with open(results + "/CFD0D/after_foaming.txt") as fl2:
            rho, dcell, xCyP0, xCO20 = loadtxt(fl2, skiprows=1, unpack=True)
        for name in names:
            xfrac[name] = 0
        xfrac['CO2'] = xCO20
        xfrac['CyP'] = xCyP0
    elif inputs['sourceOfProperty']['gasComposition'] == 'Qmom3D':
        with open(results + "/CFD3D/after_foaming.txt") as fl2:
            rho, dcell, xCyP0, xCO20 = loadtxt(fl2, skiprows=1, unpack=True)
        for name in names:
            xfrac[name] = 0
        xfrac['CO2'] = xCO20
        xfrac['CyP'] = xCyP0
    else:
        raise Exception("unknown source for gas composition")
    return xfrac


def get_fstrut(inputs, rho0):
    """Extract strut content."""
    if inputs['sourceOfProperty']['strutContent'] == 'DirectInput':
        fstrut0 = inputs['morphology']['strutContent']
    elif inputs['sourceOfProperty']['strutContent'] == 'StrutContent':
        inp = {'rho': rho0}
        model = SurrogateModel.load('strutContent')
        outputs = model.callModel(inp)
        fstrut0 = outputs['fs']
    else:
        raise Exception("unknown source for strut content")
    return fstrut0


def foamAging(flag):
    "Create initial points from foamAging.json."
    if os.getcwd().split(os.path.sep)[-1] != 'foamAging':
        raise Exception("you can initialize for foamAging only from foamAging "
                        + "folder")
    results = '../foamExpansion/results/'
    with open("./inputs/foamAging.json") as fl:
        inputs = json.load(fl)
        T0 = inputs['foamCondition']['conductivityTemperature']
        rhop = inputs['physicalProperties']['polymerDensity']
        rho0 = get_rho(inputs, results)
        dcell0 = get_dcell(inputs, results)
        xfrac = get_gas_comp(inputs, results)
        fstrut0 = get_fstrut(inputs, rho0)
    tot = 0
    for val in xfrac.itervalues():
        tot += val
    for key, val in xfrac.iteritems():
        xfrac[key] = val / tot
    eps0 = 1 - rho0 / rhop
    if flag == 1:
        ini = {
            'eps': setIP(eps0),
            'dcell': setIP(dcell0),
            'fstrut': setIP(fstrut0),
            'T': setIP(T0)
        }
        for key, val in xfrac.iteritems():
            ini["x[" + key + "]"] = setIP_fractions(val)
    elif flag == 2:
        ini = {
            'eps': setIP(eps0),
            'dcell': setIP(dcell0),
            'fstrut': setIP(fstrut0),
            'T': setIP(T0),
            'x[CO2]': [1, 0, 0, 0],
            'x[CyP]': [0, 1, 0, 0],
            'x[O2]':  [0, 0, 1, 0],
            'x[N2]':  [0, 0, 0, 1]
        }
        mat = eye(4)
        for i, key in enumerate(xfrac.iterkeys()):
            ini["x[" + key + "]"] = mat[i % 4, :].tolist()
    with open("./inputs/init_foamConductivity.json", "w") as fl:
        json.dump(ini, fl, indent=4, sort_keys=True)
    return


if __name__ == '__main__':
    if len(sys.argv) == 2:
        if sys.argv[1] == '0':
            no_refitting()
        elif sys.argv[1] == '1':
            foamAging(1)
        elif sys.argv[1] == '2':
            foamAging(2)
        else:
            my_help()
    else:
        my_help()
