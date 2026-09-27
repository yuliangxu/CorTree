# Precision-prior API

`CorTree_sampler`, `PhyloTree_sampler`, `CorTree_sampler_randominit` and
`CorTree_sampler_randominit_heldout` accept the following controls. They affect
the correlated block; independent-tree fits use `all_ind=TRUE`.

| Setting | `ghs_diag_rate` | `ghs_diag_upper` | `ghs_det_df` |
|---|---:|---:|---:|
| Flat benchmark (default) | `0` | `Inf` | `0` |
| Exponential | `log(2)` | `Inf` | `0` |
| Bounded uniform | `0` | `100` | `0` |
| Bounded determinant | `0` | `100` | fixed `a` |
| Weak soft, scale `s0` | `2*s0^2` | `Inf` | `4` |
| Strong soft, scale `s0` | `a*s0^2/2` | `Inf` | `a = N-N/K` |

Here `N` is the number of profiles and `K=n_clus` is the finite number of mixture
slots. Compute `a` once before fitting; do not update it with current occupancy.
The main report uses `s0=0.1`, `K=3` in simulation and `K=5` for DNase. These
options preserve the horseshoe off-diagonal factors and random global scale.
See the [result comparison](sampler-corrections.md) for prior definitions and
the empirical tradeoffs.

```r
N <- nrow(X)
K <- 3L
a <- N - N / K
s0 <- 0.1
fit <- CorTree::CorTree_sampler(
  X = X, n_clus = K, tree_depth = 6L, cutoff_layer = 4L,
  total_iter = 150L, burnin = 100L, init_Z = init_Z,
  c_sigma2_vec = 10, sigma_mu2 = 0.1, cov_interval = 5L,
  warm_start = 0L, ghs_det_df = a,
  ghs_diag_rate = a * s0 * s0 / 2, ghs_diag_upper = Inf
)
```

The reproduction scripts construct `X` and zero-based `init_Z`, save the
sampler-entry RNG state, and specify every study setting. The API's covariance
update default is `1`; the simulation study uses `5`, and DNase uses `3`.
The short study chains are for matched comparison, not a convergence guarantee.

For the hierarchical-scale family set `ghs_scale_hierarchy=TRUE`, with the
weak soft template above, and explicitly use the study controls
`ghs_scale_shape=3`, `ghs_scale_rate_shape=4`, `ghs_scale_rate_rate=2`.
They define `Omega[k]=t[k]*Q[k]`, `t[k]|b ~ Gamma(3,b)` and
`b ~ Gamma(4,2)` in shape/rate notation. The study controls differ from the
generic hierarchy defaults. Saved `mcmc$Sigma_inv` matrices are effective
precisions; `mcmc$ghs_scale` contains the component scales and shared rate.

Proper modes require a positive diagonal rate or a finite diagonal upper bound.
A positive determinant strength can use either support. Finite bounds and
positive rates cannot be combined in the current API. Bounded and hierarchical
fits require `warm_start=0`. The study scripts keep the ordinary assignment
updates (`z_mode=0`, `z_det_gamma0=1`).

Proper modes update empty-component precisions, disable legacy covariance
ridges/freezes, and fail on invalid numerical updates. Bounds are enforced in
the joint conditional, not by clipping draws. The determinant strength changes
the Schur-complement shape to `(n_k+a)/2+1`; for bounded priors it also changes
the integrated block mass used by the off-diagonal update. Positive diagonal
rates affect both the Schur-complement and off-diagonal conditionals.

Flat GHS remains an improper diagnostic benchmark: empty-component precision
updates are held, and legacy ridge/freeze safeguards remain. Inspect
`mcmc$ghs_prior`, `mcmc$empty_precision_updates` and
`mcmc$covariance_safeguards`. Fixed finite component truncation does not repair
an improper base measure.

## Numerical verification

Install the package into a separate library and run the actual production-kernel
checks with outputs outside the repository:

```sh
mkdir -p /path/to/R-library
R CMD INSTALL --preclean --library=/path/to/R-library .
Rscript tests/check_ghs_priors.R . /path/to/R-library /path/to/check-results
```

Install `Rcpp`, `RcppArmadillo`, `ape` and `mclust` first. The suite checks both
samplers against scalar and low-dimensional reference laws, bounded support,
strong soft settings, hierarchy updates, empty components, independent-mode
invariance and wrapper forwarding. `CORTREE_REFERENCE_LIBRARY` optionally
enables exact legacy-mode replay against a separately installed earlier build;
it is not required for standalone checks. Without arguments the script skips
these standalone checks during routine package checking. Successful numerical
tests do not establish empirical posterior convergence.
