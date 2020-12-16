    // Initialise solid field pointer lists
PtrList<coordinateSystem> coordinates(solidRegions.size());
PtrList<volScalarField> kss(solidRegions.size());
PtrList<volScalarField> fiss(solidRegions.size());
//PtrList<volScalarField> kappass(solidRegions.size());
//PtrList<fv::options> solidPotentialSources(solidRegions.size());
List<bool> residualReachedSolid(solidRegions.size(), true);
    // Populate solid field pointer lists
forAll(solidRegions, i)
{
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

	/*Info<< "    Adding fvOptions\n" << endl;
    	solidPotentialSources.set
    	(
        	i,
        	new fv::options(solidRegions[i])
    	);*/

        Info<< "    Adding to kss\n" << endl;
        kss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "ks",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

        Info<< "    Adding to fiss\n" << endl;
        fiss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fi_s",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
       );
       // just in order to check correct units
       /*Info<< "    Adding to kappass\n" << endl;
       kappass.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappas",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
		dimensionedScalar("kappas", dimensionSet(-1,-5,4,0,0,2,0), 1.0)
            )
       );*/
}

