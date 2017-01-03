"""
@file       vtkconv.py
@namespace  FoamConstruction.vtkconv
@ingroup    mod_foamConstruction
@brief      Converts binary VTK to ASCII VTK file.
@author     Pavel Ferkl
@copyright  2014-2016, MoDeNa Project. GNU Public License.
"""
import vtk
def main(filenameIn,filenameOut,origin,spacing):
    """Main function for the conversion.

    Intended for VTK files with 3D voxel data.
    Adjusts origin and spacing.
    @param[in] filenameIn binary VTK file name
    @param[in] filenameOut ASCII VTK file name
    @param[in] origin of coordinate system
    @param[in] distance between the points
    """
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
