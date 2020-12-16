fvMesh& mesh = solidRegions[i];

volScalarField& ks = kss[i];
volScalarField& fis = fiss[i];
//volScalarField& kappas = kappass[i];

//fv::options& fvOptions = solidPotentialSources[i];
bool& resReachedSolid = residualReachedSolid[i];
