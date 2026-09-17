# Technical report: corrected sampler and covariance safeguards

The sampler now uses cluster-specific sample sizes in the graphical horseshoe
update, the corrected horseshoe auxiliary and stick-breaking conditionals, and
non-truncated variance shapes. Phylogenetic updates use the actual internal-parent
counts, including for unbalanced trees. Stable binomial calculations replace
logit clipping; input validation, single-component cases, held-out draw alignment,
and post-burn-in clustering summaries are also corrected. `dahl_clustering()`
selects a sampled partition rather than averaging numeric cluster labels.

The original graphical horseshoe prior, including its flat precision-diagonal
prior, is preserved. The original numerical safeguards are retained: covariance
eigenvalues receive a small ridge when `det(precision) > 1e150`, and correlated
covariance updates stop while `det(precision) >= 1e200`. Thresholds are evaluated
on the log scale, and the internal covariance/precision pair stays synchronized.
The ridge is `1e-15` for CorTree and `10^(-11/sqrt(L))` for PhyloTree, where `L` is
the correlated block size. Inspect `fit$mcmc$covariance_safeguards` for per-component
regularization and skipped-update counts over the full chain. These guards are
numerical approximations, not exact posterior updates. Empty components retain
their correlated precision because the flat diagonal prior supplies no proper
empty-component conditional; proper mean and independent-variance priors are
still updated.

For the repository's Sim1-style experiment, explicitly set `cov_interval = 5L`.
The API default remains `1L`. A matched three-dataset diagnostic (100 samples,
K=3, depth=6, cutoff=4, 150 iterations with 100 discarded) gave corrected ARIs
0.974, 0.977, and 1.000 with interval 5, versus 0.629, 0.523, and 0.098 with interval
1. These small-run results demonstrate schedule sensitivity, not convergence or
a reproduction of the full paper study. Use multiple starts and assess mixing;
a high ARI does not rule out nearly singular covariance or frozen assignments.
The rendered example outputs in the README predate these corrections.

