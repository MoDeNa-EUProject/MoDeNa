# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:32:34 2015
Reads binary vtk file and creates ascii vtk file
@author: Pavel Ferkl
"""
import vtk
def main(filenameIn,filenameOut,origin,spacing):
    r = vtk.vtkDataSetReader()
    r.SetFileName(filenameIn)
    r.Update()
    data = vtk.vtkImageData()
    data.ShallowCopy(r.GetOutput())
    data.SetOrigin(origin)
    data.SetSpacing(spacing)
    data.Update()
    w = vtk.vtkStructuredPointsWriter()
    w.SetInput(data)
    w.SetFileName(filenameOut)
    w.Write()
