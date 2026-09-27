# Matched simulation reproduction

`simulation.R` generates the 300 Sim1 datasets and runs the selected precision
priors and four comparison methods. It needs no private files, archived results,
cluster scheduler, or `/cwork` directory. Choose an output directory outside the
source repository. The script generates no figures.

Install this checkout of CorTree and the simulation dependencies first:

```r
install.packages(c("Rcpp", "RcppArmadillo", "dplyr", "cluster", "mclust", "BiocManager"))
BiocManager::install("DirichletMultinomial")
```

```sh
mkdir -p /path/to/R-library
R CMD INSTALL --library=/path/to/R-library .
Rscript reproduction/ghs/simulation.R --output /path/to/results --library /path/to/R-library --smoke
Rscript reproduction/ghs/simulation.R --output /path/to/results --library /path/to/R-library --task 1
# Run tasks individually (including in an external array), or sequentially:
Rscript reproduction/ghs/simulation.R --output /path/to/results --library /path/to/R-library --task all
Rscript reproduction/ghs/simulation.R --output /path/to/results --summarize
```

Omit `--library` to use the active R library. A supplied library must contain the
CorTree installation actually loaded. The historical corrected runs used
R 4.4.3. Exact RNG replay also depends on the baseline package versions, especially
DirichletMultinomial; every run records versions, the CorTree binary checksum,
runner checksum, and `sessionInfo()`. Different software versions can produce
different chains even with the same seed. Do not mix versions or selected/all
settings in the same output directory: the protocol check rejects that mixture.

## Dataset and RNG protocol

Tasks 1–100 have N=200, tasks 101–200 have N=400, and tasks 201–300 have N=600.
Every task uses seed `2025 + task_id`, Mersenne-Twister/Inversion/Rejection RNGs,
1,000 bins, and two generating groups with probabilities 0.6 and 0.4. The
historical generator draws a Beta(10,10) weight and allocates exactly
`floor(total_count * weight)` observations to its first beta component. It does
not draw an independent mixture membership for each observation. Counts are
generated in the original group loop; totals are sampled from 1,000–5,000.

Initial component labels are `dplyr::ntile(rowSums(X), 3) - 1`. Kmeans on
`scale(X)` with `nstart=25`, PAM on raw counts, and a three-component
DirichletMultinomial fit run in that order **before** saving CorTree's entry RNG.
The DMM call deliberately omits its `seed` argument, matching its historical
default and RNG consumption. Flat CorTree runs next; IndTree starts from the RNG
immediately after that flat fit. Every proper-prior CorTree arm resets to the
saved flat CorTree entry RNG. Thus changing the list or order of proper arms does
not change their starting random stream or the IndTree fit.

All fitted mixtures use K=3 slots. CorTree uses depth 6, cutoff 4 (31 correlated
tree variables), `c_sigma2_vec=10`, `sigma_mu2=0.1`, no warm start, 150 iterations,
100 burn-in iterations, and covariance updates every 5 iterations. IndTree uses
cutoff 3 and `all_ind=TRUE`, with otherwise matched settings. These are short
historical chains, not a convergence recommendation.

`--smoke` uses only N=20 and 6/2 total/burn-in iterations. It writes under
`smoke/task_001` with its own protocol, and its results never enter the production
summary. It is an implementation check, not an accuracy experiment.

## Prior settings

For precision Ω, the determinant/trace factor is
`|Ω|^(a/2) exp(-rho * tr(Ω))`, combined with the graphical-horseshoe off-diagonal
prior and the positive-definite constraint. A finite upper bound additionally
requires every diagonal to be below M. All gamma specifications below use rates.

The default comparison includes four CorTree settings, plus Kmeans, PAM, DMM,
and IndTree:

| Saved method | Diagonal/determinant specification |
|---|---|
| `flat` | Original improper flat diagonal, a=0, rho=0, no upper bound; historical reference only |
| `uniform` | a=0, rho=0, M=100 |
| `soft_s0p1` | Weak determinant/trace: a=4, s0=0.1, rho=a*s0²/2=0.02 |
| `strong_soft_s0p1` | Strong determinant/trace: a=N−N/K=2N/3, s0=0.1, rho=a*s0²/2 |

Use a **different output directory** with `--all-priors` to include the other
simulation settings from the exploratory comparison:

```sh
Rscript reproduction/ghs/simulation.R --output /path/to/all-prior-results --library /path/to/R-library --task all --all-priors
```

| Additional saved method | Specification |
|---|---|
| `exponential` | a=0, rho=log(2), no upper bound |
| `bounded_4N3` | a=4N/3, rho=0, M=100 |
| `bounded_2N3` | a=2N/3, rho=0, M=100 |
| `soft_s0p3`, `soft_s1p0` | a=4, s0=0.3 or 1, rho=a*s0²/2 |
| `hier_s0p1`, `hier_s0p3`, `hier_s1p0` | Same weak prior on Q, effective precision Ω=tQ; t_k given b is Gamma(3,b), b is Gamma(4,2), with all K slots included |

The strong s0=0.3 and s0=1 arms were stopped after the initial DNase screen and
had no simulation runs in the reported study, so this runner does not include
them. The retired JMLR kappa=6 variant is also absent. Proper arms must pass
finite/SPD checks, have no covariance regularizations or skipped updates, and
emit no sampler warnings before the task is marked complete. Flat-prior
safeguards and warnings remain recorded; they do not make that prior proper.

## Artifacts, checks, and summaries

Each `simulation/task_XXX/` contains:

- `input.rds`: counts, truth, initial labels, RNG after data generation, and task metadata.
- `rng_before_cortree.rds` and `rng_before_indtree.rds`: the historical sampler entry streams.
- `Kmeans.rds`, `PAM.rds`, `DMM.rds`: baseline fits and warnings.
- A directory for each Bayesian method with the fit (including precision and
  any hierarchy traces), entry/exit RNG, Dahl assignments, occupancy and
  precision diagnostics, mixing CSV, and warning log.
- `metrics.csv`, `metadata.rds`, `prior_settings.csv`, and `sessionInfo.txt`.
- `completion.rds` with artifact MD5 checksums and a final `COMPLETE` marker
  bound to that manifest. MD5 here detects accidental corruption; it is not a
  security signature. Failed attempts retain `FAILED.txt` and partial artifacts.

The primary ARI uses the Dahl partition over all 50 retained allocations. Flat,
uniform, exponential, bounded, and IndTree results retain the historical package
matrix implementation. Weak, hierarchical, and strong soft-prior results retain
the later exact integer contingency implementation. Both minimize the same Dahl
loss; their tie decisions can differ at floating-point precision. Metadata names
the implementation. Integer label averages are not a clustering estimator and
are not used as the primary result.

The per-fit diagnostics include canonicalized distinct partitions (invariant to
component relabeling), adjacent partition agreement, occupied-slot counts,
positive definiteness, diagonal extrema, and covariance safeguards. The saved
log-likelihood is a conditional kernel, not a predictive score or marginal
likelihood. An unchanged retained partition does not establish convergence.

`--summarize` needs only base R. It verifies completion markers and artifact
checksums, then writes `summary/task_status.csv`, `summary/metrics.csv`,
`summary/simulation_summary.csv`, and `summary/STATUS.txt`. A final mean and SD
for a given N require **all 100** matched replicates. Partial means appear only
in explicitly named `partial_*` columns, with status `INCOMPLETE`. A paired
difference against the flat rerun and its Monte Carlo standard error are
reported only for a complete 100-replicate group. Smoke output is excluded.

Already complete tasks are verified and skipped. An incomplete task is rerun
deterministically and its own partial files are replaced. A `.running` lock
prevents two processes from running the same task simultaneously; after a killed
process, remove its lock only once you have confirmed that it stopped. Concurrent
tasks may share the output root; run a final `--summarize` after they finish.
Full completion demonstrates that the configured computation finished and
passed its numerical checks. It does not make these exploratory, retrospectively
selected priors an independent validation study.
