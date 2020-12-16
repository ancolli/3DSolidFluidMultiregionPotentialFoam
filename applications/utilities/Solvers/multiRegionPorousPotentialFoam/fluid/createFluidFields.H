// Initialise fluid field pointer lists
PtrList<volScalarField> kfs(fluidRegions.size());
PtrList<volScalarField> fifs(fluidRegions.size());
//PtrList<volScalarField> kappafs(fluidRegions.size());
//PtrList<fv::options> fluidFvOptions(fluidRegions.size());
List<bool> residualReachedFluid(fluidRegions.size(), true);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to kfs\n" << endl;
        kfs.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kf",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to fifs\n" << endl;
        fifs.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fi_f",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
       );

       // just in order to check correct units
       /*Info<< "    Adding to kappafs\n" << endl;
       kappafs.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappaf",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
		dimensionedScalar("kappaf", dimensionSet(-1,-5,4,0,0,2,0), 1.0)
            )
       );*/

       /*Info<< "    Adding fvOptions\n" << endl;
       fluidFvOptions.set
       (
           i,
           new fv::options(fluidRegions[i])
       );*/
}
