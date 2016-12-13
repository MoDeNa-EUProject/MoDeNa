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
