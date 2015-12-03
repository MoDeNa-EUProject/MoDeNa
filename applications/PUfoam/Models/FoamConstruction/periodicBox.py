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
Normals=[
    [[-1,0,0],[0,-1,0],[0,0,-1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,-1,0],[0,0,-1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,-1,0],[0,0,-1],[+1,0,0],[0,-1,0],[0,0,-1]],
    [[-1,0,0],[0,+1,0],[0,0,-1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,-1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,-1],[+1,0,0],[0,-1,0],[0,0,-1]],
    [[-1,0,0],[0,+1,0],[0,0,-1],[-1,0,0],[0,+1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,-1],[-1,0,0],[0,+1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,-1],[+1,0,0],[0,+1,0],[0,0,-1]],

    [[-1,0,0],[0,-1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,-1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,-1,0],[0,0,+1],[+1,0,0],[0,-1,0],[0,0,-1]],
    [[-1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[+1,0,0],[0,-1,0],[0,0,-1]],
    [[-1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,+1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,+1,0],[0,0,-1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[+1,0,0],[0,+1,0],[0,0,-1]],

    [[-1,0,0],[0,-1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,+1]],
    [[+1,0,0],[0,-1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,+1]],
    [[+1,0,0],[0,-1,0],[0,0,+1],[+1,0,0],[0,-1,0],[0,0,+1]],
    [[-1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,+1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,-1,0],[0,0,+1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[+1,0,0],[0,-1,0],[0,0,+1]],
    [[-1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,+1,0],[0,0,+1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[-1,0,0],[0,+1,0],[0,0,+1]],
    [[+1,0,0],[0,+1,0],[0,0,+1],[+1,0,0],[0,+1,0],[0,0,+1]],
        ]
regions=[]
for j in xrange(2):
    polydata=reader
    for i in xrange(6):
        plane=vtk.vtkPlane()
        plane.SetOrigin(Origins[i])
        plane.SetNormal(Normals[j][i])
        # Clip it
        clipper = vtk.vtkClipPolyData()
        clipper.SetInputConnection(polydata.GetOutputPort())
        clipper.SetClipFunction(plane)
        polydata=clipper
        polydata.Update()
    regions.append(vtk.vtkPolyData())
    regions[j].ShallowCopy(polydata.GetOutput())
# Move it
transform = vtk.vtkTransform()
transform.Translate([1,0,0])
transformFilter = vtk.vtkTransformPolyDataFilter()
transformFilter.SetTransform(transform)
transformFilter.SetInputConnection(polydata.GetOutputPort())
transformFilter.Update()
regions[1].ShallowCopy(transformFilter.GetOutput())
# Append the two meshes
appendFilter = vtk.vtkAppendPolyData()
if vtk.VTK_MAJOR_VERSION <= 5:
    for j in xrange(2):
        appendFilter.AddInputConnection(regions[j].GetProducerPort())
else:
    for j in xrange(2):
        appendFilter.AddInputData(regions[j])
appendFilter.Update()
#  Remove any duplicate points.
cleanFilter = vtk.vtkCleanPolyData()
cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
cleanFilter.Update()
# Final data to be saved and displayed
finalData=cleanFilter
# Write the stl file to disk
stlWriter = vtk.vtkSTLWriter()
stlWriter.SetFileName(filenameClipped)
stlWriter.SetInputConnection(finalData.GetOutputPort())
stlWriter.Write()
# Create mappper and actor for rendering
mapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(finalData.GetOutput())
else:
    mapper.SetInputConnection(finalData.GetOutputPort())
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
