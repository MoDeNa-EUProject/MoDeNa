"""
@file
Converts binary VTK to ASCII VTK file.

@author    Pavel Ferkl
@copyright 2014-2016, MoDeNa Project. GNU Public License.
@ingroup   foam_constr
"""
## Main function for the conversion.
#
#  Intended for VTK files with 3D voxel data.
#  Adjusts origin and spacing.
import vtk
def main(filenameIn,filenameOut,origin,spacing):
    r = vtk.vtkDataSetReader()
    r.SetFileName(filenameIn)
    r.Update()
    data = vtk.vtkImageData()
    data.ShallowCopy(r.GetOutput())
    data.SetOrigin(origin)
    data.SetSpacing(spacing)
    w = vtk.vtkStructuredPointsWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        w.SetInputConnection(data.GetProducerPort())
    else:
        w.SetInputData(data)
    w.SetFileName(filenameOut)
    w.Write()
