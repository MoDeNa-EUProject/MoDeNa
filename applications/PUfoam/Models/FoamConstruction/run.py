"""
@author: Pavel Ferkl
"""
import os
import FoamGeometryConstruction_Periodic
import periodicBox
import vtkconv

vtkconv.main("PeriodicRVEBox.vtk","PeriodicRVEBox-ascii.vtk")

#
# ########## Create periodic RVE
# # Foam parameters
# MU=0.15
# SIGMA=0.01
# NumOfCells=8
# # Set jobs
# packing=1
# tesselation=1
# geometry=1
# statistics=0
# hypermesh=0
# deleteFiles=1
# filename="PeriodicRVE"
# FoamGeometryConstruction_Periodic.main(MU,SIGMA,NumOfCells,filename,packing,\
#     tesselation,geometry,statistics,hypermesh,deleteFiles)
# ########## Convert .geo to .stl
# os.system("gmsh -2 -format stl "+filename+".geo")
# # Move foam to periodic box
# filenameIn = filename+".stl"
# filename = filename+"Box"
# filenameOut = filename+".stl"
# # Define domain origin and size
# xmin=4
# ymin=4
# zmin=4
# dx=4
# dy=4
# dz=4
# render=0
# periodicBox.main(filenameIn,filenameOut,xmin,ymin,zmin,dx,dy,dz,render)
# ########## Convert .stl to .ply
# os.system("meshconv "+filename+".stl -c ply")
# ########## Binarize and save as .vtk
# if os.path.isfile(filename+'.vtk'):
#     os.remove(filename+'.vtk')
# os.system("binvox -e -rotz -rotx -rotz -rotz -t vtk "+filename+".ply")
# ########## Convert binary .vtk to ascii .vtk
# filenameIn = filename+".vtk"
# filename = filename+"-ascii"
# filenameOut = filename+".vtk"
# vtkconv.main(filenameIn,filenameOut)
