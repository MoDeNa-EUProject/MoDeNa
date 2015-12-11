# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:32:34 2015
Reads binary vtk file and creates ascii vtk file
@author: Pavel Ferkl
"""
from __future__ import division
import vtk
def main(filenameIn,filenameOut,dx,dy,dz,vx,vy,vz):
    r = vtk.vtkDataSetReader()
    r.SetFileName(filenameIn)
    r.Update()
    data = vtk.vtkImageData()
    data.ShallowCopy(r.GetOutput())
    data.SetOrigin(dx,dy,dz)
    data.SetSpacing(dx/vx,dy/vy,dz/vz)
    data.Update()
    #w = vtk.vtkDataSetWriter()
    w = vtk.vtkStructuredPointsWriter()
    # w.SetInputConnection(data.GetProducerPort())
    w.SetInput(data)
    w.SetFileName(filenameOut)
    w.Write()
