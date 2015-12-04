# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:32:34 2015
Reads binary vtk file and creates ascii vtk file
@author: Pavel Ferkl
"""
import vtk
def main(filenameIn,filenameOut):
    r = vtk.vtkDataSetReader()
    r.SetFileName(filenameIn)
    r.Update()
    data = vtk.vtkImageData()
    data.ShallowCopy(r.GetOutput())
    data.SetOrigin(4,4,4)
    data.SetSpacing(0.015625,0.015625,0.015625)
    data.Update()
    #w = vtk.vtkDataSetWriter()
    w = vtk.vtkStructuredPointsWriter()
    # w.SetInputConnection(data.GetProducerPort())
    w.SetInput(data)
    w.SetFileName(filenameOut)
    w.Write()
