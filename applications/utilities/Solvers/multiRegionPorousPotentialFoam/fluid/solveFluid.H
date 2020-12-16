{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix fifEqn
        (
               fvm::laplacian(kf, fif) 
           /*)
	  ==
    		fvOptions(kappaf, fif)*/        
        );

        fifEqn.relax();
	//fvOptions.constrain(fifEqn);
        eqnResidual = fifEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
	//fvOptions.correct(fif);    
    }
}


Info<< "Min/max fif:" << min(fif).value() << ' '
    << max(fif).value() << endl;
