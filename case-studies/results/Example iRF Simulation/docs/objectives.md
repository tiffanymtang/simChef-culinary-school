This example asks whether iterative random forests (iRF) can recover a known, high-order interaction when the covariates and response noise resemble a real genomics application. It is a compact adaptation of the enhancer simulations in Section S5.4 of [Basu et al. (2018)](https://doi.org/10.1073/pnas.1711236115).

The experiment has three objectives:

1. Measure held-out classification performance as the training sample size increases

2. Measure recovery of known interactions

3. Compare the stability assigned to true and false interactions
