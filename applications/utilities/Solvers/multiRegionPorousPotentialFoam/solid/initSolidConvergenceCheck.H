// initialize values for convergence checks

scalar eqnResidual = 1, maxResidual = 0;
scalar convergenceCriterion_fi_s = 0;

simple.readIfPresent("conv_fi_s", convergenceCriterion_fi_s);
