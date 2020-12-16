// check convergence
resReachedPorous = false;
if (maxResidual_fi_f < convergenceCriterion_fi_f && maxResidual_fi_s < convergenceCriterion_fi_s && maxResidual_C_ < convergenceCriterion_C_)
{
    Info<< "reached convergence criterion: fi_f in Porous: " << convergenceCriterion_fi_f << endl;
    Info<< "reached convergence criterion: fi_s in Porous: " << convergenceCriterion_fi_s << endl;
    Info<< "reached convergence criterion: C_ in Porous: " << convergenceCriterion_C_ << endl;

    resReachedPorous = true;
}
