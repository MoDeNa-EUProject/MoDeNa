/**
* @ingroup mod_3Dcfd
* @file checkGellingPoint.H
* @brief checks if the gelling point reaches.
*
* @detail
* This file counts the number of cells filled with foam and the
* number of cells when the gelling criterion has been satisfied.
* If the ratio is greater than 0.5, it is assumed that the gelling points
* is reached.
*/
///@cond
scalar foamCellCounter = 0.0;
scalar gelCellCounter = 0.0;

forAll(mesh.C(), celli)
{
    if (alpha2[celli] > 0.5)
    {
        foamCellCounter += 1;
    }
    if (XOH[celli] >= XOH_Gel)
    {
        gelCellCounter += 1;
    }
}
if ((gelCellCounter/foamCellCounter) > 0.5)
{
    gellingPoint = true;
}
else
{
    gellingPoint = false;
}
///@endcond