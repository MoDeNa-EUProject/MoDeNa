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