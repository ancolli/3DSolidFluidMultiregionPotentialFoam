    
   // Added support for cases without porous regions
    const wordList porousNames
    (
        rp.found("porous") ? rp["porous"] : wordList(0)
    );


   // const wordList porousNames(rp["porous"]);

    PtrList<fvMesh> porousRegions(porousNames.size());

    forAll(porousNames, i)
    {
        Info<< "Create porous mesh for region " << porousNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        porousRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    porousNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
