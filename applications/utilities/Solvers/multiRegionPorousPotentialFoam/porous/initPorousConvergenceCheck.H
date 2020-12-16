// initialize values for convergence checks

scalar eqnResidual_fi_s = 1, maxResidual_fi_s = 0;
scalar eqnResidual_fi_f = 1, maxResidual_fi_f = 0;
scalar eqnResidual_C_ = 1, maxResidual_C_ = 0;

scalar convergenceCriterion_fi_s = 0;
scalar convergenceCriterion_fi_f = 0;
scalar convergenceCriterion_C_ = 0;

simple.readIfPresent("conv_fi_s", convergenceCriterion_fi_s);
simple.readIfPresent("conv_fi_f", convergenceCriterion_fi_f);
simple.readIfPresent("conv_C_", convergenceCriterion_C_);
