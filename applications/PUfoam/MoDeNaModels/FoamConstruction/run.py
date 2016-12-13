#!/usr/bin/env python
"""
@file
Python script, which organizes creation of the foam.

@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   foam_constr
"""
from __future__ import division
import os
import sys
from blessings import Terminal
import json
import datetime
from scipy.optimize import minimize_scalar as minimize_scalar
import FoamGeometryConstruction_Periodic
import periodicBox
import vtkconv
dx=dy=dz=1 # size of RVE
dedge=6
########## Read input file
with open('input.json') as data_file:
    data = json.load(data_file)
locals().update(data) # Creates variables from dictionary
filenameBox=filename+"Box"
filenameBoxAscii = filenameBox+"-ascii"
filenameBoxStruts = filenameBox+"Struts"
filenameDescriptors = "descriptors.out"
filenameParameters = "parameters.out"
########## Create terminal for colour output
term = Terminal()
## Objective function for finding size of box, which would give desired porosity
def porOpt(vx):
    vx=int(vx)
    vy=vx
    vz=vx
    if os.path.isfile(filenameBox+'.vtk'):
        os.remove(filenameBox+'.vtk')
    if not os.path.isfile(filenameBox+'.ply'):
        raise SystemError(".ply file is missing. Nothing to binarize.")
    os.system("binvox -e -d {0:d} -rotz -rotx -rotz -rotz -t vtk ".format(vx)+
        filenameBox+".ply >binvox.out")
    with open('binvox.out') as data_file:
        for line in data_file:
            if "counted" in line:
                solidVoxel,totalVoxel=\
                    [int(s) for s in line.split() if s.isdigit()]
                eps=1-solidVoxel/totalVoxel
                print "dimension: {0:4d}, porosity: {1:f}".format(vx,eps)
                return (eps-porosity)**2
## Objective function for finding size of box, which
#
#  would give desired porosity and strut content.
def porfsOpt(x):
    global dedge
    vx=int(x)
    vy=vx
    vz=vx
    if os.path.isfile(filenameBox+'.vtk'):
        os.remove(filenameBox+'.vtk')
    if not os.path.isfile(filenameBox+'.ply'):
        raise SystemError(".ply file is missing. Nothing to binarize.")
    os.system("binvox -e -d {0:d} -rotz -rotx -rotz -rotz -t vtk ".format(vx)+
        filenameBox+".ply >binvox.out")
    filenameIn = filenameBox+".vtk"
    filenameOut = filenameBoxAscii+".vtk"
    origin=[dx,dy,dz]
    spacing=[dx/vx,dy/vy,dz/vz]
    vtkconv.main(filenameIn,filenameOut,origin,spacing)
    f=open("foamreconstr.in","w")
    f.write("0\n")
    f.write("1\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("0\n")
    f.write("{0:f}\n".format(dedge))
    f.write("{0:f}\n".format(1-strutContent*(1-porosity)))
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
    f.write(filenameBoxStruts+"\n")
    f.write(filenameBoxAscii+".vtk\n")
    f.write(filename+".gnu\n")
    f.write("name\n")
    f.write(filenameDescriptors+"\n")
    f.write(filenameParameters+"\n")
    f.close()
    os.system("./foamreconstr/foamreconstr")
    f=open(filenameDescriptors,"r")
    eps=float(f.readline())
    fs=float(f.readline())
    f.close()
    f=open(filenameParameters,"r")
    dedge=float(f.readline())
    f.close()
    resid=((eps-porosity)/porosity)**2
    print "dimension: {0:4d}, porosity: {1:f}".format(vx,eps)+\
        ", strut content: {0:f}".format(fs)
    return resid
ts = datetime.datetime.now()
########## Create periodic RVE of foam
print(
    term.yellow +
    "Create periodic RVE of foam" +
    term.normal
)
FoamGeometryConstruction_Periodic.main(MU,SIGMA,NumOfCells,filename,packing,
    alternativePackingAlgorithm,tesselation,visualizeTesselation,geometry,
    statistics,hypermesh,deleteFiles,dx,dy,dz)
if moveToPeriodicBox:
    ########## Convert .geo to .stl
    print(
        term.yellow +
        "Convert .geo to .stl" +
        term.normal
    )
    os.system("gmsh -n -2 -format stl "+filename+".geo >gmsh.out")
    if deleteFiles:
        os.remove("gmsh.out")
    ########## Move to periodic box
    print(
        term.yellow +
        "Move to periodic box" +
        term.normal
    )
    filenameIn = filename+".stl"
    filenameOut = filenameBox+".stl"
    xmin=dx
    ymin=dy
    zmin=dz
    periodicBox.main(filenameIn,filenameOut,xmin,ymin,zmin,dx,dy,dz,renderBox)
    if deleteFiles:
        os.remove(filenameIn)
    ########## Convert .stl to .ply
    print(
        term.yellow +
        "Convert .stl to .ply" +
        term.normal
    )
    os.system("meshconv "+filenameBox+".stl -c ply")
    if deleteFiles:
        os.remove(filenameBox+'.stl')
if binarizeBox:
    ########## Binarize and save as .vtk
    if strutContent==0:
        # Find the size of box, which would give desired porosity
        # This method is not optimal, since the solver doesn't know that the
        # function takes only integer arguments
        print(
            term.yellow +
            "Optimizing porosity" +
            term.normal
        )
        res=minimize_scalar(porOpt,bracket=[100,120],method='Brent',tol=1e-2)
        vx=res.x
        vx=int(vx)
        print 'box size: {0:d}'.format(vx)
        vy=vx
        vz=vx
        print(
            term.yellow +
            "Creating and saving optimal foam" +
            term.normal
        )
        porOpt(vx) # Call it with the optimized box size
        if deleteFiles:
            os.remove("binvox.out")
            os.remove(filename+".ply")
        ######### Convert binary .vtk to ascii .vtk
        print(
            term.yellow +
            "Convert binary .vtk to ascii .vtk" +
            term.normal
        )
        filenameIn = filenameBox+".vtk"
        filenameOut = filenameBoxAscii+".vtk"
        origin=[dx,dy,dz]
        spacing=[dx/vx,dy/vy,dz/vz]
        vtkconv.main(filenameIn,filenameOut,origin,spacing)
        if deleteFiles:
            os.remove(filenameIn)
    else:
        print(
            term.yellow +
            "Optimizing porosity and strut content" +
            term.normal
        )
        res=minimize_scalar(porfsOpt,bracket=[150,200],method='Brent',tol=1e-2)
        # res=minimize_scalar(porfsOpt,bounds=[200,250],method='bounded',tol=2e0)
        vx=res.x
        vx=int(vx)
        print 'optimal box size: {0:d}'.format(vx)
        vy=vx
        vz=vx
        print(
            term.yellow +
            "Creating and saving optimal foam" +
            term.normal
        )
        if os.path.isfile(filenameBox+'.vtk'):
            os.remove(filenameBox+'.vtk')
        if not os.path.isfile(filenameBox+'.ply'):
            raise SystemError(".ply file is missing. Nothing to binarize.")
        os.system("binvox -e -d {0:d}".format(vx)+" -rotz -rotx -rotz -rotz "+
            "-t vtk "+filenameBox+".ply >binvox.out")
        filenameIn = filenameBox+".vtk"
        filenameOut = filenameBoxAscii+".vtk"
        origin=[dx,dy,dz]
        spacing=[dx/vx,dy/vy,dz/vz]
        vtkconv.main(filenameIn,filenameOut,origin,spacing)
        f=open("foamreconstr.in","w")
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write("0\n")
        f.write("1\n")
        f.write("0\n")
        f.write("{0:f}\n".format(dedge))
        f.write("{0:f}\n".format(1-strutContent*(1-porosity)))
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
        f.write(filenameBoxStruts+"\n")
        f.write(filenameBoxAscii+".vtk\n")
        f.write(filename+".gnu\n")
        f.write("name\n")
        f.write(filenameDescriptors+"\n")
        f.write(filenameParameters+"\n")
        f.close()
        os.system("./foamreconstr/foamreconstr")
        if deleteFiles:
            os.remove(filenameBoxAscii+".vtk")
            os.remove(filenameDescriptors)
            os.remove(filenameParameters)
            os.remove("binvox.out")
            os.remove(filenameBox+".vtk")
            os.remove("foamreconstr.in")
tf = datetime.datetime.now()
te = tf - ts
print "Foam created in: ",te
