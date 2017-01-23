#!/usr/bin/env python
import sys
import os
import json
from numpy import loadtxt
from modena import SurrogateModel
def my_help():
    "print usage of the module"
    print('Usage: prep_init_foamConductivity v')
    print('v: integer')
    print("0 - try to initialize with points so that foamConductivity "+
        "does not go out of range")
    print('1 - initialize with points from foamAging.json')
    print('2 - initialize with points from foamAging.json except for gas '
        +'composition - foamConductivity should not go out of range in '
        +'foamAging simulation')
    return
def no_refitting():
    """
    initialize with dummy data that should ensure that foamConductivity
    does not go out of range
    """
    ini={
        "eps": [0.9,0.0,0.96,0.99,0.7,0.5],
        "dcell": [200e-6,0.0,300e-6,100e-6,1e-2,200e-6],
        "fstrut": [0.0,1.0,0.7,0.6,0.0,0.9],
        "T": [280,549,300,350,330,300],
        "x[CO2]": [1,0,0,0,1,0],
        "x[CyP]": [0,1,0,0,0,1],
        "x[O2]": [0,0,1,0,0,0],
        "x[N2]": [0,0,0,1,0,0]
    }
    with open("./inputs/init_foamConductivity.json","w") as fl:
        json.dump(ini,fl,indent=4,sort_keys=True)
    return
def setIP(a0):
    "creates a list with small perturbation around the input value"
    a=[]
    for i in xrange(4):
        a.append(a0)
    upar=1-1e-4
    opar=1+1e-4
    a[0]=a[0]*upar
    if a0==0:
        a[1]=1e-4
    else:
        a[1]=a[1]*opar
    return a
def setIP_fractions(a0):
    "creates a list with small perturbation around the input value"
    a=[]
    for i in xrange(4):
        a.append(a0)
    upar=1-1e-4
    opar=1+1e-4
    if a0>1:
        a0=1.0
    a[0]=a[0]*upar
    if a0==0:
        a[1]=1e-4
    else:
        a[1]=a[1]*opar
        if a[1]>1:
            a[1]=1.0
    return a
def get_rho(inputs, results):
    if inputs['sourceOfProperty']['foamDensity']=='DirectInput':
        rho0=inputs['morphology']['foamDensity']
    elif inputs['sourceOfProperty']['foamDensity']=='BubbleGrowth':
        with open(results+"/bubbleGrowth/after_foaming.txt") as fl2:
            rho0,dcell,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
    elif inputs['sourceOfProperty']['foamDensity']=='Qmom0D':
        with open(results+"/CFD0D/after_foaming.txt") as fl2:
            rho0,dcell,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
    elif inputs['sourceOfProperty']['foamDensity']=='Qmom3D':
        with open(results+"/CFD3D/after_foaming.txt") as fl2:
            rho0,dcell,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
    else:
        raise Exception("unknown source for foam density")
    return rho0
def get_dcell(inputs, results):
    if inputs['sourceOfProperty']['cellSize']=='DirectInput':
        dcell0=inputs['morphology']['cellSize']
    elif inputs['sourceOfProperty']['cellSize']=='BubbleGrowth':
        with open(results+"/bubbleGrowth/after_foaming.txt") as fl2:
            rho,dcell0,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
    elif inputs['sourceOfProperty']['cellSize']=='Qmom0D':
        with open(results+"/CFD0D/after_foaming.txt") as fl2:
            rho,dcell0,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
    elif inputs['sourceOfProperty']['cellSize']=='Qmom3D':
        with open(results+"/CFD3D/after_foaming.txt") as fl2:
            rho,dcell0,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
    else:
        raise Exception("unknown source for cell size")
    return dcell0
def get_gas_comp(inputs, results):
    if inputs['sourceOfProperty']['gasComposition']=='DirectInput':
        xO20=inputs['foamCondition']['initialComposition']['O2']
        xN20=inputs['foamCondition']['initialComposition']['N2']
        xCO20=inputs['foamCondition']['initialComposition']['CO2']
        xCyP0=inputs['foamCondition']['initialComposition']['Cyclopentane']
    elif inputs['sourceOfProperty']['gasComposition']=='BubbleGrowth':
        with open(results+"/bubbleGrowth/after_foaming.txt") as fl2:
            rho,dcell,xCyP0,xCO20=loadtxt(fl2,skiprows=1,unpack=True)
            xO20=xN20=0
    elif inputs['sourceOfProperty']['gasComposition']=='Qmom0D':
        with open(results+"/CFD0D/after_foaming.txt") as fl2:
            rho,dcell,xCyP0,xCO20=loadtxt(fl2,skiprows=1,unpack=True)
            xO20=xN20=0
    elif inputs['sourceOfProperty']['gasComposition']=='Qmom3D':
        with open(results+"/CFD3D/after_foaming.txt") as fl2:
            rho,dcell,xCyP0,xCO20=loadtxt(fl2,skiprows=1,unpack=True)
            xO20=xN20=0
    else:
        raise Exception("unknown source for gas composition")
    return xO20,xN20,xCO20,xCyP0
def get_fstrut(inputs, results):
    if inputs['sourceOfProperty']['strutContent']=='DirectInput':
        fstrut0=inputs['morphology']['strutContent']
    elif inputs['sourceOfProperty']['strutContent']=='StrutContent':
        inputs = {'rho': rho0}
        model=SurrogateModel.load('strutContent')
        outputs=model.callModel(inputs)
        fstrut0=outputs['fs']
    else:
        raise Exception("unknown source for strut content")
    return fstrut0
def foamAging(flag):
    "create initial points from foamAging.json"
    if os.getcwd().split(os.path.sep)[-1] != 'foamAging':
        raise Exception("you can initialize for foamAging only from foamAging "
            +"folder")
    results='../foamExpansion/results/'
    with open("./inputs/foamAging.json") as fl:
        inputs=json.load(fl)
        T0=inputs['foamCondition']['conductivityTemperature']
        rhop=inputs['physicalProperties']['polymerDensity']
        rho0=get_rho(inputs, results)
        dcell0=get_dcell(inputs, results)
        xO20,xN20,xCO20,xCyP0=get_gas_comp(inputs, results)
        fstrut0=get_fstrut(inputs, results)
    s=xO20+xN20+xCO20+xCyP0
    xO20=xO20/s
    xN20=xN20/s
    xCO20=xCO20/s
    xCyP0=xCyP0/s
    eps0=1-rho0/rhop
    if flag==1:
        ini = {
            'eps': setIP(eps0),
            'dcell': setIP(dcell0),
            'fstrut': setIP(fstrut0),
            'T': setIP(T0),
            'x[CO2]': setIP_fractions(xCO20),
            'x[CyP]': setIP_fractions(xCyP0),
            'x[O2]': setIP_fractions(xO20),
            'x[N2]': setIP_fractions(xN20),
        }
    elif flag==2:
        ini = {
            'eps': setIP(eps0),
            'dcell': setIP(dcell0),
            'fstrut': setIP(fstrut0),
            'T': setIP(T0),
            'x[CO2]': [1,0,0,0],
            'x[CyP]': [0,1,0,0],
            'x[O2]':  [0,0,1,0],
            'x[N2]':  [0,0,0,1]
        }
    with open("./inputs/init_foamConductivity.json","w") as fl:
        json.dump(ini,fl,indent=4,sort_keys=True)
    return
if __name__=='__main__':
    n=len(sys.argv)
    if n==2:
        if sys.argv[1]=='0':
            no_refitting()
        elif sys.argv[1]=='1':
            foamAging(1)
        elif sys.argv[1]=='2':
            foamAging(2)
        else:
            my_help()
    else:
        my_help()
