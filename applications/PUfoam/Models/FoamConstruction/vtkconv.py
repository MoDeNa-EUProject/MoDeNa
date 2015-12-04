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

    #w = vtk.vtkDataSetWriter()
    w = vtk.vtkStructuredPointsWriter()
    w.SetInput( r.GetOutput() )
    w.SetFileName(filenameOut)
    w.Write()
