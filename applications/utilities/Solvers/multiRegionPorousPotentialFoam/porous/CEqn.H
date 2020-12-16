{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
	    
        //forAll(speciesNameList, i)
        //{

	    dimensionedScalar a("a",dimensionSet(0,2,-1,0,0,0,0), 7.35e-10);//[kg m s K mol A cd]
	    dimensionedScalar b("b",dimensionSet(0,3.5,-1,0,-0.5,0,0), -1.676e-11);//[kg m s K mol A cd]
	    D_ = a+b*pow(C_,0.5);

	    fvScalarMatrix CEqn
            (
                fvm::ddt(C_)
              + fvm::div(phi, C_)
              - fvm::laplacian(D_, C_)
             ==
                fvOptions(C_)
            );

            CEqn.relax();
            fvOptions.constrain(CEqn);

            eqnResidual_C_ = CEqn.solve().initialResidual();
	    maxResidual_C_ = max(eqnResidual_C_, maxResidual_C_);

	    fvOptions.correct(C_);

	    // Bound
	    C_.max(SMALL);
	    C_.correctBoundaryConditions();    
	//}
    }
}

Info<< "Min/max C_:" << min(C_).value() << ' '
    << max(C_).value() << endl;




// ************************************************************************* //
