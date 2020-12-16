{
//dimensionedScalar cte("cte",dimensionSet(0,-3,0,0,0,1,0), 2000.0);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
	fisporous.storePrevIter();
        fvScalarMatrix fisporousEqn//solve//
        (
                 fvm::laplacian(ksporous, fisporous) 
		==
		 fvOptions(kappap,fisporous)//fvOptions(kappas, fisporous) //-cte//
        );
	//fisporous.relax();
        fisporousEqn.relax();

        fvOptions.constrain(fisporousEqn);
	
	eqnResidual_fi_s = fisporousEqn.solve().initialResidual();
        maxResidual_fi_s = max(eqnResidual_fi_s, maxResidual_fi_s); 

	fvOptions.correct(fisporous);

	Jvol == fvc::laplacian(ksporous, fisporous);
    }

}

Info<< "Min/max fisporous:" << min(fisporous).value() << ' '
    << max(fisporous).value() << endl;



