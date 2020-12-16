
int nFluidControlled = fluidRegions.size();
int nPorousControlled = porousRegions.size();
int nSolidControlled = solidRegions.size();


    int nFluidConv = 0;
    int nPorousConv = 0;
    int nSolidConv = 0;

    // Sum of all converged regions (Note: if no residual control is used
    // the residualReached* flag is already set to true)
    forAll(fluidRegions, i)
    {
        if (residualReachedFluid[i])
        {
            nFluidConv++;
        }
    }

    forAll(porousRegions, i)
    {
        if (residualReachedPorous[i])
        {
            nPorousConv++;
        }
    }


    forAll(solidRegions, i)
    {
        if (residualReachedSolid[i])
        {
            nSolidConv++;
        }
    }

    if (nFluidConv == nFluidControlled && nPorousConv == nPorousControlled && nSolidConv == nSolidControlled)
    {
        // Activate flag to go to the 'Final' loop using the 'Final'
        // relaxation factors
        allRegionsConverged = true;
    }


if (!runTime.run())
{
    Info<< "\nRegions not converged after " << runTime.timeName()
        << " iterations" << endl;
}
else if (runTime.run() && allRegionsConverged)
{
    Info<< "\nRegions converged after " << runTime.timeName()
        << " iterations" << endl;

    // Save converged time step and end the run
    runTime.writeAndEnd();
}
