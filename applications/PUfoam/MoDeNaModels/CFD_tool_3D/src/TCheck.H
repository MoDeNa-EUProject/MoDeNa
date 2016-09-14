forAll(mesh.C(), celli)
{
    if (T[celli] > 500.0)
    {
        T[celli] = 500.0;
    }
}
