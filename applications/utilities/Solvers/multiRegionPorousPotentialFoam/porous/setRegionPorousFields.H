fvMesh& mesh = porousRegions[i];

volScalarField& ksporous = ksporouss[i];
volScalarField& kfporous = kfporouss[i];
volScalarField& fisporous = fisporouss[i];
volScalarField& fifporous = fifporouss[i];
volScalarField& C_ = C_s[i];
volScalarField& D_ = D_s[i];
//volVectorField& U = Us[i]; //not necessary 
surfaceScalarField& phi = phiFluid[i];

volScalarField& Jvol = Jvols[i];

volScalarField& kappap = kappaps[i];

fv::options& fvOptions = porousFvOptions[i];

bool& resReachedPorous = residualReachedPorous[i];

