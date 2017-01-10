/**
* @ingroup mod_3Dcfd
* @file conversionCheck.H
* @brief checks the boundedness of coversion values (for polyol and water).
*/
///@cond
forAll(mesh.C(), celli)
{
    if (XW[celli] < 0.0 || XOH[celli] > XOH_Gel || XW[celli] > 1.0)
    {
        XW[celli] = ROOTVSMALL;
    }
    if (XOH[celli] < 0.0 || XOH[celli] > 1.0)
    {
        XOH[celli] = ROOTVSMALL;
    }
    if (alpha2[celli] < 0.5)
    {
        XOH[celli] = ROOTVSMALL;
        XW[celli] = ROOTVSMALL;
    }
    if (XOH[celli] > XOH_Gel)
    {
        XOH[celli] = XOH_Gel;
    }
}
///@endcond