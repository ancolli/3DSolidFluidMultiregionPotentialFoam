    // Added support for cases without solid regions
    const wordList solidNames
    (
        rp.found("solid") ? rp["solid"] : wordList(0)
    );

    PtrList<fvMesh> solidRegions(solidNames.size());

    forAll(solidNames, i)
    {
        Info<< "Create solid mesh for region " << solidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        solidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    solidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

