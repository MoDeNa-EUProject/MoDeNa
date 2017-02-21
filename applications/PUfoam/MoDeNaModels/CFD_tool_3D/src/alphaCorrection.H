/**
@ingroup mod_3Dcfd
@file alphaCorrection.H
@brief checks the boundedness of alpha (phase volume fraction)

* @param gasvalue volume fraction of gas phase
* @param foamvalue volume fraction of foam phase
* @param alphasum summation of volume fraction for two phases
* @param corrector constant value
* @param foamold volume fraction of foam phase at old time
* @param gasold volume fraction of gas phase at old time
*/
double gasvalue, foamvalue,alphasum,corrector,foamold,gasold;
corrector = 0.0001;
forAll(mesh.C(), celli)
{

    gasvalue    = alpha1[celli];
    foamvalue   = alpha2[celli];

    foamold     = alpha2.oldTime()[celli];
    gasold      = alpha1.oldTime()[celli];
    alphasum    = gasvalue + foamvalue;

    while (alphasum > 1.0)
    {
        if (gasvalue > foamvalue)
        {
            gasvalue = gasvalue - (Foam::mag(gasvalue - gasold))*corrector;
        }
        else
        {
            foamvalue = foamvalue - (Foam::mag(foamvalue - foamold))*corrector;
        }
    }
}
