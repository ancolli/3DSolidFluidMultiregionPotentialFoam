// initialize values for convergence checks

scalar eqnResidual = 1, maxResidual = 0;
scalar convergenceCriterion_fi_f = 0;

simple.readIfPresent("conv_fi_f", convergenceCriterion_fi_f);
