The covariates come from `case-studies/data/enhancer.Rdata`, which contains 7,987 observations and 80 enhancer measurements. The simulation retains columns 46 through 80, corresponding to 35 transcription-factor features. Related measurements are assigned to the same feature group using `varnames.all`; for example, `kr1` and `kr2` belong to the `kr` group.

The noiseless response is the locally spiky sparse Boolean rule

$$
f(X)=I\!\left[
(kr1>1.25\ \text{or}\ kr2>1.25)
\land(twi1>1.25\ \text{or}\ twi2>1.25)
\land(D1>1.25)
\land(hb1>1.25\ \text{or}\ hb2>1.25)
\land(wt\_ZLD>75)
\right].
$$

Thus the true grouped interaction is `D1_hb_kr_twi_wt.ZLD`. Its subsets define the true lower-order interactions used by the recovery evaluator.

RF-derived noise is estimated once before running the simulation. For each of 20 repetitions, the original enhancer observations are split 50/50 within response class. A classification forest is fit to one half using `randomForest::randomForest()` defaults, and class-1 probabilities are predicted only for the held-out half. For observation $i$, $p_i$ is the mean of all held-out predictions available for that observation.

Each simulation replicate draws a new response

$$
Y_i=I\{f(X_i)=1\ \text{or}\ U_i<p_i\}, \qquad U_i\sim\text{Uniform}(0,1).
$$

This construction preserves the known LSS signal while adding heterogeneous positive-class structure learned from the real data. Finally, observations are sampled within response class: approximately `n` observations form the training set and all remaining observations form the test set. The DGP returns the training and test matrices, factor responses with levels 0 and 1, feature-group labels, and the five true feature groups.
