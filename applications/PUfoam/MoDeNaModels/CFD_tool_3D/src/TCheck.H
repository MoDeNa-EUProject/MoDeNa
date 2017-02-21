/**
* @ingroup mod_3Dcfd
* @file TCheck.H
* @brief checks the boundedness of temperature values
* 
* @detail
* This file checks the boundedness of calculated temperature in the foam phase.
* The criterion is selected based on physical evidences.
*/
forAll(mesh.C(), celli)
{
    if (T[celli] > 500.0)
    {
        T[celli] = 500.0;
    }
}
