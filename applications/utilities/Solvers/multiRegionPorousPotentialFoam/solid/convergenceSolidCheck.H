// check convergence
resReachedSolid = false;
if (maxResidual < convergenceCriterion_fi_s)
{
    Info<< "reached convergence criterion: fi_s in Solid: " << convergenceCriterion_fi_s << endl;

    resReachedSolid = true;
}
