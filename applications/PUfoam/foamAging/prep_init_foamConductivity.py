#!/usr/bin/env python
import sys
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
def foamAging():
    "create initial points from foamAging.json"
    results='../foamExpansion/results/'
    with open("./inputs/foamAging.json") as fl:
        inputs=json.load(fl)
        T0=inputs['foamCondition']['conductivityTemperature']
        rhop=inputs['physicalProperties']['polymerDensity']
        if inputs['sourceOfProperty']['foamDensity']=='DirectInput':
            rho0=inputs['morphology']['foamDensity']
        elif inputs['sourceOfProperty']['foamDensity']=='BubbleGrowth':
            with open(results+"/bubbleGrowth/after_foaming.txt") as fl2:
                rho0,dcell,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
        elif inputs['sourceOfProperty']['foamDensity']=='Qmom0D':
            with open(results+"/CFD0D/after_foaming.txt") as fl2:
                rho0,dcell,var,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
        elif inputs['sourceOfProperty']['foamDensity']=='Qmom3D':
            with open(results+"/CFD3D/after_foaming.txt") as fl2:
                rho0,dcell,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
        else:
            raise Exception("unknown source for foam density")
        if inputs['sourceOfProperty']['cellSize']=='DirectInput':
            dcell0=inputs['morphology']['cellSize']
        elif inputs['sourceOfProperty']['cellSize']=='BubbleGrowth':
            with open(results+"/bubbleGrowth/after_foaming.txt") as fl2:
                rho,dcell0,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
        elif inputs['sourceOfProperty']['cellSize']=='Qmom0D':
            with open(results+"/CFD0D/after_foaming.txt") as fl2:
                rho,dcell0,var,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
        elif inputs['sourceOfProperty']['cellSize']=='Qmom3D':
            with open(results+"/CFD3D/after_foaming.txt") as fl2:
                rho,dcell0,xCyP,xCO2=loadtxt(fl2,skiprows=1,unpack=True)
        else:
            raise Exception("unknown source for cell size")
        if inputs['sourceOfProperty']['gasComposition']=='DirectInput':
            xAir0=inputs['foamCondition']['initialComposition']['Air']
            xCO20=inputs['foamCondition']['initialComposition']['CO2']
            xCyP0=inputs['foamCondition']['initialComposition']['Cyclopentane']
        elif inputs['sourceOfProperty']['gasComposition']=='BubbleGrowth':
            with open(results+"/bubbleGrowth/after_foaming.txt") as fl2:
                rho,dcell,xCyP0,xCO20=loadtxt(fl2,skiprows=1,unpack=True)
                xAir0=0
        elif inputs['sourceOfProperty']['gasComposition']=='Qmom0D':
            with open(results+"/CFD0D/after_foaming.txt") as fl2:
                rho,dcell,var,xCyP0,xCO20=loadtxt(fl2,skiprows=1,unpack=True)
                xAir0=0
        elif inputs['sourceOfProperty']['gasComposition']=='Qmom3D':
            with open(results+"/CFD3D/after_foaming.txt") as fl2:
                rho,dcell,xCyP0,xCO20=loadtxt(fl2,skiprows=1,unpack=True)
                xAir0=0
        else:
            raise Exception("unknown source for gas composition")
        if inputs['sourceOfProperty']['strutContent']=='DirectInput':
            fstrut0=inputs['morphology']['strutContent']
        elif inputs['sourceOfProperty']['strutContent']=='StrutContent':
            inputs = {'rho': rho0}
            model=SurrogateModel.load('strutContent')
            outputs=model.callModel(inputs)
            fstrut0=outputs['fs']
        else:
            raise Exception("unknown source for strut content")
    eps0=1-rho0/rhop
    ini = {
        'eps': setIP(eps0),
        'dcell': setIP(dcell0),
        'fstrut': setIP(fstrut0),
        'T': setIP(T0),
        'x[CO2]': setIP(xCO20),
        'x[CyP]': setIP(xCyP0),
        'x[O2]': setIP(xAir0*0.21),
        'x[N2]': setIP(xAir0*0.79),
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
            foamAging()
        else:
            my_help()
    else:
        my_help()
