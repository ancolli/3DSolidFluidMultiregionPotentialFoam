    // Initialise porous field pointer lists

PtrList<volScalarField> ksporouss(porousRegions.size());
PtrList<volScalarField> kfporouss(porousRegions.size());
PtrList<volScalarField> fisporouss(porousRegions.size());
PtrList<volScalarField> fifporouss(porousRegions.size());
PtrList<volScalarField> Jvols(porousRegions.size());
PtrList<volScalarField> C_s(porousRegions.size());
PtrList<volScalarField> D_s(porousRegions.size());
PtrList<volVectorField> Us(porousRegions.size());
PtrList<surfaceScalarField> phiFluid(porousRegions.size());

PtrList<volScalarField> kappaps(porousRegions.size());

PtrList<fv::options> porousFvOptions(porousRegions.size());

List<bool> residualReachedPorous(porousRegions.size(), true);

    // Populate 3D field pointer lists
forAll(porousRegions, i)
{
        Info<< "*** Reading Porous mesh thermophysical properties for region "
            << porousRegions[i].name() << nl << endl;


        Info<< "    Adding to ksporouss\n" << endl;
        ksporouss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "ks",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousRegions[i]
            )
        );

	Info<< "    Adding to kfporouss\n" << endl;
        kfporouss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kf",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousRegions[i]
            )
        );


        Info<< "    Adding to fisporouss\n" << endl;
        fisporouss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fi_s",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousRegions[i]
            )
        );

	Info<< "    Adding to fifporouss\n" << endl;
        fifporouss.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "fi_f",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousRegions[i]
            )
       );

       //creates object to be able to have the J*Ae output in paraview
	Info<< "    Adding J*Ae\n" << endl;

        Jvols.set
        (
            i,
            new volScalarField
            (
                 IOobject
                 (
                     "Jvol",
                     runTime.timeName(),
                     porousRegions[i],
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 porousRegions[i],
                 dimensionedScalar("Jvol", dimensionSet(0,-3,0,0,0,1,0), 0.0)//[kg m s K mol A cd]
             )
	);


       Info<< "    Adding to C_s\n" << endl;
        C_s.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "C_",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousRegions[i]
            )
       );

       Info<< "    Adding to D_s\n" << endl;
        D_s.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "D_",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                porousRegions[i]
            )
       );

       Info<< "Reading field U\n" << endl;

    	Us.set
        (
            i,
            new volVectorField 
    	    (
        	IOobject
        	(
            	    "U",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE//AUTO_WRITE
        	),
                porousRegions[i]
            )
	);

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                porousRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE//AUTO_WRITE
            ),
            fvc::flux(Us[i])
        )
    );

       // just in order to check correct units
       Info<< "    Adding Artificial units\n" << endl;
       kappaps.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kappap",
                    runTime.timeName(),
                    porousRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                porousRegions[i],
		dimensionedScalar("kappap", dimensionSet(-1,-5,4,0,0,2,0), 1.0)
            )
       );


       //fvOptions

       Info<< "    Adding fvOptions\n" << endl;
       porousFvOptions.set
       (
           i,
	   new fv::options(porousRegions[i])

       );

}

