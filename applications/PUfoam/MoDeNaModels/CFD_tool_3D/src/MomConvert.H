/**
* @ingroup mod_3Dcfd
* @file MomConvert.H
* @brief converts moments based on the unit volume of foam
*/

///@cond
if ( nNodes == 2 )
{
    M0 = mZero*(1.0 - (mOne/(1.0+mOne)));
    M1 = mOne*(1.0 - (mOne/(1.0+mOne)));
    M2 = mTwo*(1.0 - (mOne/(1.0+mOne)));
    M3 = mThree*(1.0 - (mOne/(1.0+mOne)));
}
if ( nNodes == 3 )
{
    M0 = mZero*(1.0 - (mOne/(1.0+mOne)));
    M1 = mOne*(1.0 - (mOne/(1.0+mOne)));
    M2 = mTwo*(1.0 - (mOne/(1.0+mOne)));
    M3 = mThree*(1.0 - (mOne/(1.0+mOne)));
    M4 = mFour*(1.0 - (mOne/(1.0+mOne)));
    M5 = mFive*(1.0 - (mOne/(1.0+mOne)));
}
///@endcond