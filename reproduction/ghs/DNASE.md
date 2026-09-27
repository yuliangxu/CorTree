# Reproducing the DNase comparisons

`dnase.R` runs the corrected CorTree and IndTree comparisons for one dataset.
It does not reproduce the historical executable or promise identical random
draws across R/compiler versions. The archived studies used R 4.4.3.
Install the current package first, then run the script from a checkout.
All output must be outside that checkout; no input data or figures are included.

```sh
R CMD INSTALL --library=/path/to/R-library .
Rscript reproduction/ghs/dnase.R --data-dir /path/to/DNA_data \
  --library /path/to/R-library --output /path/to/results --dataset REST \
  --init rowsum --prior soft_s0p1 --smoke
```

The library directory must already exist. Omit `--library` to use the ordinary
R library search path. `--help` prints options without fitting.
`--out` is an alias for `--output`.

## External inputs

Supply `--data-dir`, or set `CORTREE_DNA_DATA`. For each requested dataset,
the directory must contain these two RDS files (replace `REST` with `NRF1`):

- `REST.K562.DNase.counts.mat.rds`: nonnegative integer count matrix. Its first
  half of columns contains one strand and its second half the matching other
  strand, with aligned rows and genomic bins.
- `REST.K562.sites.chip.labels.tss.dist.rds`: a data frame with one row per count
  matrix row and columns `pwm.score`, `strand`, `tss.dist`, `chip_label`, `chip`.
  `chip_label` is the binary ChIP reference used only for evaluation.

The runner adds corresponding strand columns, then keeps positive-strand sites
with total count at least 50 and PWM score at least 13, preserving input order.
It requires the report's dimensions: REST 4,551 × 220; NRF1 8,148 × 211.
Use the processed data described in the
[paper's data-availability section](https://arxiv.org/html/2509.15480v2).
Raw ENCODE downloads alone do not satisfy this input contract.

## Selected report settings

```sh
Rscript reproduction/ghs/dnase.R --data-dir /path/to/DNA_data \
  --library /path/to/R-library --output /path/to/results --dataset REST
Rscript reproduction/ghs/dnase.R --data-dir /path/to/DNA_data \
  --library /path/to/R-library --output /path/to/results --dataset NRF1
```

Each command defaults to `--prior selected --init all`: corrected flat,
bounded uniform, weak soft at 0.1, strong soft at 0.1, and IndTree, each with
count-sum, PWM and TSS-distance starts (15 fits per dataset). Choose one start
with `--init rowsum`, `pwm`, or `tss_dist`, and one prior by name to limit work.

| Prior name | Precision-prior settings |
|---|---|
| `flat` | Flat precision diagonals; improper benchmark with legacy safeguards |
| `uniform` | Diagonals in `(0,100)`; random global horseshoe scale |
| `soft_s0p1` | `a=4`, `rho=0.02`, no upper bound |
| `strong_soft_s0p1` | Fixed `a=N-N/5`, `rho=0.005*a`, no upper bound |
| `indtree` | Independent tree: `all_ind=TRUE`, no correlated precision prior |

All tree fits use five slots, depth 9, correlated cutoff 3 for REST or 4 for
NRF1, `c_sigma2_vec=10`, `sigma_mu2=0.1`, no warm start, and 150 iterations
with 100 discarded. CorTree updates precision every three iterations; IndTree
updates independent variances every iteration. A label-invariant Dahl partition
summarizes all 50 retained allocation draws.

Initialization uses the package's deterministic quantile bins, including its
handling of tied quantiles. The RNG is reset to Mersenne-Twister/Inversion/Rejection
and seed 2025 immediately before **each** sampler call. No random operation is
inserted between this reset and fitting; the entry/exit RNG states are saved.
This matches the archived sampler-entry protocol independently of prior/start
execution order. Exact reproduction additionally requires the same processed
input rows and compatible R, package, compiler and numerical libraries.

`--smoke` computes initial labels on the full filtered data, retains the first
64 rows, and uses 6 iterations with 2 discarded. Strong prior strength still
uses full filtered N. Smoke results occupy a separate output directory and
must not be compared with the report's production ARIs.

## Other tested settings

`--prior all` (or `--all-priors`) explicitly requests all listed settings for the chosen dataset;
it can be expensive (48 REST or 42 NRF1 fits with all three starts). In addition
to the selected settings, its menu contains:

- `exponential`: `a=0`, `rho=log(2)`, infinite upper bound.
- `bounded_r0.25`, `bounded_r1`, `bounded_r4`: upper bound 100 and fixed
  `a=r*N/5`. The `r=0` control is already represented by `uniform`.
- `soft_s0p3`, `soft_s1p0`: `a=4`, `rho=2*s0^2`.
- `hier_s0p1`, `hier_s0p3`, `hier_s1p0`: the preceding soft template,
  `Omega=t*Q`, `t|b ~ Gamma(3,b)`, `b ~ Gamma(4,2)` (shape/rate).
- REST only: `strong_soft_s0p3`, `strong_soft_s1p0`, with fixed `a=N-N/5`
  and `rho=a*s0^2/2`. These scales were not run on NRF1 in the archived study.

This standalone runner does **not** implement the earlier sequential screening
controller. Running NRF1 is an explicit independent request; it neither checks
nor asserts that REST passed a performance gate. The original strong-soft study
continued only scale 0.1 after its REST screen. These are exploratory benchmarks,
and data reuse, invariant chains and initialization sensitivity limit inference.

## Non-tree baselines

Install `cluster`, `CENTIPEDE`, and Bioconductor's `DirichletMultinomial` separately,
then request only the four non-tree methods:

```sh
Rscript reproduction/ghs/dnase.R --data-dir /path/to/DNA_data \
  --library /path/to/R-library --output /path/to/results --dataset REST --prior baselines
```

K-means uses `scale(X)`, two clusters and 25 starts; PAM uses raw `X` and two
clusters. CENTIPEDE uses the intercept and PWM score, thresholding posterior
binding probability at 0.5. DMM uses two components and an explicit seed 2025.
Each method resets its own RNG, independently of tree fits. `--init` does not
apply to these methods. These reproduce the corrected baselines: notably, REST
DMM's corrected ARI differs from the paper's unseeded historical baseline.

## Outputs and reruns

Each fit writes `metrics.csv`, `assignments.csv`, `fit.rds`, metadata, saved
initial labels/RNG states, numerical diagnostics, warnings and session information
under `OUT/production/DATASET/PRIOR_INIT` (or `OUT/smoke/...`). File checksums and
a `COMPLETE` marker are written only after fitting and validation finish.
Each invocation also writes a `summary_PRIOR_INIT.csv` table for its requested fits.
Reissuing the same command verifies completed outputs and skips them; a changed
input/package/configuration or incomplete directory stops rather than overwrites
it. Inspect failures and use a fresh output root for a retry. Warnings are retained
in the metrics and `warnings.txt`; inspect them before interpreting results.
Numerical validity and a stable retained partition do not establish convergence.
