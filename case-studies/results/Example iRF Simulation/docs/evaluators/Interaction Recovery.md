The known active feature groups are $S=\{D1,hb,kr,twi,wt.ZLD\}$. For interaction order $s\in\{2,3,4,5\}$, the evaluator enumerates all $\binom{5}{s}$ subsets of $S$. A target subset is counted as recovered when every target feature occurs in at least one interaction returned by iRF; the returned interaction may contain additional features.

The order-specific recovery rate is

$$
\text{Recovery}(s)=
\frac{\text{number of true order-}s\text{ interactions recovered}}
{\binom{5}{s}}.
$$

Consequently, the denominators for orders 2 through 5 are 10, 10, 5, and 1. If iRF returns no interactions, recovery is defined as zero at every order. simChef averages the replicate-level recovery rates within each sample-size and interaction-order combination.
