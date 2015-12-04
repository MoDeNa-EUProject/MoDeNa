# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 12:25:27 2015
Takes the representative volume element and moves it into a box with
periodic boundary conditions
@author: Pavel Ferkl
"""
import vtk
# print vtk.VTK_MAJOR_VERSION # Check the version
# Define input and output files
filename = "PeriodicRVE.stl"
filenameOut = "PeriodicRVEDomain.stl"
# Define domain size
dx=4
dy=4
dz=4
# Read the file and create polydata
reader = vtk.vtkSTLReader()
reader.SetFileName(filename)
# Define planes for clipping
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
# Define directions for moving clipped regions
Direction=[
    [dx,dy,dz],
    [0,dy,dz],
    [-dx,dy,dz],
    [dx,0,dz],
    [0,0,dz],
    [-dx,0,dz],
    [dx,-dy,dz],
    [0,-dy,dz],
    [-dx,-dy,dz],
    [dx,dy,0],
    [0,dy,0],
    [-dx,dy,0],
    [dx,0,0],
    [0,0,0],
    [-dx,0,0],
    [dx,-dy,0],
    [0,-dy,0],
    [-dx,-dy,0],
    [dx,dy,-dz],
    [0,dy,-dz],
    [-dx,dy,-dz],
    [dx,0,-dz],
    [0,0,-dz],
    [-dx,0,-dz],
    [dx,-dy,-dz],
    [0,-dy,-dz],
    [-dx,-dy,-dz],
]
regions=[]
n=27
for j in xrange(n):
    polydata=reader
    # Clip it with all 6 planes
    for i in xrange(6):
        plane=vtk.vtkPlane()
        plane.SetOrigin(Origins[i])
        plane.SetNormal(Normals[j][i])
        clipper = vtk.vtkClipPolyData()
        clipper.SetInputConnection(polydata.GetOutputPort())
        clipper.SetClipFunction(plane)
        polydata=clipper
        polydata.Update()
    # Move it if not empty
    if polydata.GetOutput().GetLength()>0:
        transform = vtk.vtkTransform()
        transform.Translate(Direction[j])
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(polydata.GetOutputPort())
        transformFilter.Update()
        regions.append(vtk.vtkPolyData())
        regions[j].ShallowCopy(transformFilter.GetOutput())
    else:
        regions.append(vtk.vtkPolyData())
        regions[j].ShallowCopy(polydata.GetOutput())
# Append the all regions
appendFilter = vtk.vtkAppendPolyData()
if vtk.VTK_MAJOR_VERSION <= 5:
    for j in xrange(n):
        appendFilter.AddInputConnection(regions[j].GetProducerPort())
else:
    for j in xrange(n):
        appendFilter.AddInputData(regions[j])
appendFilter.Update()
#  Remove any duplicate points
cleanFilter = vtk.vtkCleanPolyData()
cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
cleanFilter.Update()
# Final data to be saved and displayed
finalData=cleanFilter
# Write the stl file to disk
stlWriter = vtk.vtkSTLWriter()
stlWriter.SetFileName(filenameOut)
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
