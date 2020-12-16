{
//dimensionedScalar cte("cte",dimensionSet(0,-3,0,0,0,1,0), 2000.0);
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {

        fifporous.storePrevIter();

        fvScalarMatrix fifporousEqn//solve //
        (
                 fvm::laplacian(kfporous, fifporous)
		==
		 fvOptions(kappap,fifporous)//  cte//    
        );

	fifporous.relax();

        fifporousEqn.relax();

	fvOptions.constrain(fifporousEqn);
	
	eqnResidual_fi_f = fifporousEqn.solve().initialResidual();
        maxResidual_fi_f = max(eqnResidual_fi_f, maxResidual_fi_f);

	fvOptions.correct(fifporous);
    }
}


Info<< "Min/max fifporous:" << min(fifporous).value() << ' '
    << max(fifporous).value() << endl;
