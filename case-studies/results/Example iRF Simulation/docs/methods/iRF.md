The method fits [iterative random forests](https://doi.org/10.1073/pnas.1711236115) with `iRF::iRF()`. iRF repeatedly refits a random forest, using feature importance from the preceding iteration to update the feature-sampling probabilities. Candidate interactions are extracted from forest decision paths with random intersection trees and assigned bootstrap stability scores.

Five iterations are fit (`niter = 5`); unsigned interactions are used (`signed = FALSE`); computation uses one core; and 20 bootstrap samples are used to estimate stability.

The random-intersection-tree configuration is `depth = 5`, `ntree = 100`, `nchild = 2`, `class.id = 1`, and `min.nd = 1`. Other iRF and random-forest arguments use the installed package defaults. For each iteration, the method returns the RF fit, discovered interactions, held-out response, and known true feature groups needed by the evaluators.
