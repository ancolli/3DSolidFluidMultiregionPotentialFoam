{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix fisEqn
        (
            //(
                fvm::laplacian(ks, fis)
		//==
    		   //fvOptions(kappas, fis)
	    //)       
        );
        fisEqn.relax();
	//fvOptions.constrain(fisEqn);	
	eqnResidual = fisEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
	//fvOptions.correct(fis); 
    }
}

Info<< "Min/max fis:" << min(fis).value() << ' '
    << max(fis).value() << endl;
