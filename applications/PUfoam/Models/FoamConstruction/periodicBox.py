# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 12:25:27 2015

@author: Pavel Ferkl
"""
import vtk
# print vtk.VTK_MAJOR_VERSION # Check the version
# Define input and output files
filename = "PeriodicRVE.stl"
filenameClipped = "PeriodicRVEClipped.stl"
# Read the file and create polydata
reader = vtk.vtkSTLReader()
reader.SetFileName(filename)
# Create plane for clipping
Origins=[[4,4,4],[4,4,4],[4,4,4],[8,8,8],[8,8,8],[8,8,8]]
Normals=[[1,0,0],[0,1,0],[0,0,1],[-1,0,0],[0,-1,0],[0,0,-1]]
planes=[]
for i in xrange(6):
    planes.append(vtk.vtkPlane())
    planes[i].SetOrigin(Origins[i])
    planes[i].SetNormal(Normals[i])
# Clip it
clipper = vtk.vtkClipPolyData()
clipper.SetInputConnection(reader.GetOutputPort())
clipper.SetClipFunction(planes[0])
clipper.GenerateClipScalarsOn()
clipper.GenerateClippedOutputOn()
clipper.SetValue(0.5)
# Write the stl file to disk
stlWriter = vtk.vtkSTLWriter()
stlWriter.SetFileName(filenameClipped)
stlWriter.SetInputConnection(clipper.GetOutputPort())
stlWriter.Write()
# Create mappper and actor for rendering
mapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(clipper.GetOutput())
else:
    mapper.SetInputConnection(clipper.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
# Create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
# Create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
# Assign actor to the renderer
ren.AddActor(actor)
# Enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()
