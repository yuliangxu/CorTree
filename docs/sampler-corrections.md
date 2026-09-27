# CorTree precision priors: simulation and DNase-seq comparison

September 27, 2026

**Proper priors can preserve CorTree's strong simulation performance, but none of the selected settings recovers the paper's performance across all DNase initializations.** Bounded uniform and weak soft GHS perform well in simulation. Strong soft GHS improves DNase clustering, especially REST with count-sum initialization, at a cost in simulation accuracy.

This report compares five CorTree versions and the paper's competing methods. **Paper** means the published result, not a new run of the historical executable; **rerun** means a completed fit using the current comparison protocol. The paper's baseline is **K-means**, rather than KNN. PAM is partitioning around medoids; DMM is a Dirichlet multinomial mixture. Published values come from [Tables 1–2 of the CorTree paper, arXiv v2](https://arxiv.org/html/2509.15480v2).

## 1. Prior settings

The package implements the proper alternatives reported here while retaining flat GHS as the default for compatibility. See the [prior API](precision-prior-api.md) for explicit settings and the [reproduction guide](../reproduction/ghs/README.md) for portable simulation and DNase scripts. Published paper rows remain historical references, not outputs of the corrected sampler.

Let $\Omega=\Sigma^{-1}$ be a component's latent-logit precision matrix. Each prior has joint kernel $\mathbf{1}_{\{\Omega\succ0\}}w(\Omega)H(\Omega,\lambda,\tau)$, where $H$ contains zero-mean Gaussian off-diagonal factors with variance $\lambda_{ij}^2\tau^2$ and standard half-Cauchy local/global scales. The table specifies the remaining factor.

| CorTree version | Factor $w(\Omega)$ and setting | Proper prior? |
|---|---|---|
| Paper version | Flat diagonal factor $1$; historical sampler | No |
| Corrected flat GHS | Same factor $1$; corrected sampler | No |
| Bounded uniform GHS | $\prod_j\mathbf{1}_{\{0<\omega_{jj}<100\}}$ | Yes |
| Weak soft GHS, $s_0=0.1$ | $\det(\Omega)^{a/2}\exp\{-\rho\operatorname{tr}(\Omega)\}$; $a=4$, $\rho=0.02$ | Yes |
| Strong soft GHS, $s_0=0.1$ | Same soft factor; $a=N-N/K$, $\rho=0.005a$ | Yes |

Both soft settings are unbounded and satisfy $E(\Sigma_{jj})=s_0^2=0.01$; $s_0$ is the square root of the prior mean **marginal variance**. Strong soft changes prior concentration as well as correlation regularization. Its strength is fixed throughout each fit: $a=2N/3$ for simulation ($K=3$), and $a=4N/5$ for DNase ($K=5$). This is a sample-size-dependent prior schedule.

The uniform bound is part of the prior support, not clipping of unconstrained draws. Proper variants update empty-component precisions and disable the legacy ridge/freeze rules. Flat GHS has no normalizable empty-component precision conditional; finite mixture truncation does not repair this. It is retained here as a diagnostic benchmark.

## 2. Simulation

There are **100 matched datasets at each $N=200,400,600$**, with two generating groups and three fitted components. Profiles contain 1,000 bins. Tree fits use depth 6, 150 iterations, 100 burn-in, and count-sum initialization; CorTree has 31 correlated splits and updates covariance every five iterations. The four CorTree reruns reuse identical datasets and initial labels. Primary rerun partitions use Dahl's label-invariant summary.

Entries are **mean adjusted Rand index (across-dataset SD)**; higher ARI is better. Published values retain their original precision; reruns are rounded to three decimals.

| Method / setting | Source | $N=200$ | $N=400$ | $N=600$ |
|---|---|---|---|---|
| CorTree, paper version | Paper | 0.94 (0.13) | 0.96 (0.10) | 0.96 (0.13) |
| CorTree, corrected flat | Rerun | 0.969 (0.057) | 0.957 (0.152) | 0.978 (0.104) |
| CorTree, bounded uniform | Rerun | 0.967 (0.103) | 0.971 (0.109) | 0.981 (0.102) |
| CorTree, weak soft $s_0=0.1$ | Rerun | 0.969 (0.068) | 0.962 (0.127) | 0.977 (0.103) |
| CorTree, strong soft $s_0=0.1$ | Rerun | 0.891 (0.110) | 0.885 (0.126) | 0.938 (0.088) |
| IndTree | Paper | 0.89 (0.06) | 0.91 (0.04) | 0.91 (0.04) |
| IndTree, corrected | Rerun | 0.839 (0.075) | 0.878 (0.060) | 0.883 (0.050) |
| K-means | Paper | 0.30 (0.04) | 0.30 (0.03) | 0.30 (0.02) |
| K-means | Rerun | 0.302 (0.037) | 0.300 (0.030) | 0.303 (0.023) |
| PAM | Paper | 0.34 (0.05) | 0.34 (0.04) | 0.35 (0.03) |
| PAM | Rerun | 0.342 (0.051) | 0.345 (0.044) | 0.349 (0.034) |
| DMM | Paper | 0.81 (0.15) | 0.76 (0.15) | 0.76 (0.15) |
| DMM | Rerun | 0.810 (0.154) | 0.755 (0.155) | 0.759 (0.151) |

**Finding.** Uniform and weak soft retain the high mean ARI of corrected flat GHS, without an improper precision prior. Strong soft has lower means at every sample size. Its mean exceeds corrected IndTree at each size, but the $N=400$ difference is only 0.0069 (paired dataset MCSE 0.0148), so superiority there is not established. K-means, PAM and DMM reruns remain close to their published simulation results.

## 3. DNase-seq: REST and NRF1

The main datasets contain **4,551 REST** and **8,148 NRF1** sites, filtered to PWM score $\ge13$, positive motif strand and total count $\ge50$. Tree fits use five component slots, depth 9, and 150 iterations with 100 burn-in; CorTree updates covariance every three iterations. Count, PWM and TSS refer to initialization by count sum, motif score and distance to the nearest transcription start site. Each cell is one fit, not a repeated-seed mean. ARI compares all-site partitions against ChIP binding labels, an external biological proxy; no post-hoc cluster merging is used.

| Tree method / setting | REST count | REST PWM | REST TSS | NRF1 count | NRF1 PWM | NRF1 TSS |
|---|---:|---:|---:|---:|---:|---:|
| CorTree, paper version | 0.319 | 0.418 | 0.244 | 0.199 | 0.203 | 0.240 |
| CorTree, corrected flat | 0.0574 | 0.0724 | 0.0482 | 0.1168 | 0.1684 | 0.1742 |
| CorTree, bounded uniform | 0.0576 | 0.0681 | 0.0494 | 0.1161 | 0.1657 | 0.1846 |
| CorTree, weak soft $s_0=0.1$ | 0.0570 | 0.0702 | 0.0473 | 0.1167 | 0.1679 | 0.1770 |
| CorTree, strong soft $s_0=0.1$ | 0.3349 | 0.1966 | 0.1067 | 0.1968 | 0.1934 | 0.2162 |
| IndTree, paper | 0.054 | 0.075 | 0.043 | 0.115 | 0.157 | 0.176 |
| IndTree, corrected | 0.0569 | 0.0706 | 0.0419 | 0.1135 | 0.1586 | 0.1785 |

The other methods fit two clusters and do not share these three tree initializations. CENTIPEDE uses PWM as a model covariate.

| Baseline | REST paper | REST rerun | NRF1 paper | NRF1 rerun |
|---|---:|---:|---:|---:|
| K-means | 0.069 | 0.0687 | 0.005 | 0.0052 |
| PAM | 0.217 | 0.2168 | 0.061 | 0.0611 |
| CENTIPEDE (PWM) | 0.336 | 0.3364 | 0.089 | 0.0887 |
| DMM | 0.012 | 0.1061 | 0.000 | -0.0001 |

The REST DMM discrepancy is preserved explicitly. DMM is unaffected by CorTree's sampler corrections; the rerun fixed its seed to 2025, whereas the historical baseline was not explicitly seeded. Its change cannot be attributed to the CorTree fixes.

**REST finding.** Strong soft improves all three initializations over corrected flat, uniform and weak soft. Its count-sum ARI, 0.3349, exceeds the paper's corresponding 0.319 and is close to CENTIPEDE's 0.3364. However, its PWM/TSS values remain substantially below the paper's 0.418/0.244; both also fall below PAM and CENTIPEDE. Uniform and weak soft do not recover the paper's REST advantage.

**NRF1 finding.** Strong soft is the best of the selected corrected CorTree settings for every initialization. Its ARIs are within 0.024 of the corresponding paper values and exceed corrected IndTree and all four non-tree baselines. These are observed point estimates, not evidence of statistical equivalence to the paper.

## 4. Interpretation and limits

- **The historical gap is not solely a prior effect.** Corrected flat retains the same improper diagonal factor yet loses much of the paper's DNase performance. The old GHS update used total $N$ with component-specific scatter; the correction uses component size $n_k$. Strong soft is motivated by the associated precision inflation, but is a coherent new prior, not a reconstruction of the historical sampler. Other sampler corrections and partition-summary changes also matter; attribution to one change would require an ablation.
- **IndTree needs no additional rerun for these tables.** The corrected 300 simulations and six main DNase fits already exist. IndTree skips GHS, but shared stick-breaking, independent-variance, tree-indexing, latent-logit and empty-component corrections affected it. Its published and corrected results are therefore shown separately.
- **Convergence remains unresolved.** All 24 selected CorTree DNase fits and all six corrected IndTree fits retain a single partition across their 50 saved draws, despite differing across starts. Numerical validity and high ARI do not establish convergence. Simulation chains also frequently have invariant partitions. Longer diagnostics for earlier settings did not resolve DNase initialization dependence; strong soft has not received that longer-run assessment.
- **These are exploratory comparisons.** The selected reruns comprise 1,200 completed CorTree simulation fits and 24 DNase fits reused from completed studies. Prior development and REST-first screening used these same benchmarks; passing the 0.1-ARI continuation tolerance does not establish equivalence. The matched simulation generator follows the archived code's fixed within-profile component counts, which differ from the paper's description of independent mixture draws. Historical paper rows are reference values, not paired reruns.

**Recommendation.** Use weak soft $s_0=0.1$ as the unbounded proper alternative that preserves simulation performance, with bounded uniform as a simple proper comparator. Present strong soft $s_0=0.1$ as the promising DNase compromise, while retaining its simulation cost and REST initialization sensitivity. The current evidence supports this tradeoff, not one uniformly superior prior.

This report focuses on simulation and DNase. American Gut comparisons cover only corrected flat and uniform, have no external clustering truth, and cannot support a five-version accuracy comparison.

The results were consolidated from the completed corrected-reproduction, bounded-uniform, weak-soft and strong-soft studies. Per-fit artifacts and publication figures are retained separately from this public repository. The [sampler correction commit](https://github.com/yuliangxu/CorTree/commit/8b92cf3) records the changes underlying the corrected benchmark.

The rendered example outputs in the README predate those sampler corrections and are not part of the results above.
