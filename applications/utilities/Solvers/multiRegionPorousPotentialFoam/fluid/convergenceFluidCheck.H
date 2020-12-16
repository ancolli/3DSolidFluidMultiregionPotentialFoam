// check convergence
resReachedFluid = false;
if (maxResidual < convergenceCriterion_fi_f)
{
    Info<< "reached convergence criterion: fi_f in fluid: " << convergenceCriterion_fi_f << endl;
    
    resReachedFluid = true;
}
