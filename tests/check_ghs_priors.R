#!/usr/bin/env Rscript
# Standalone production-kernel verification; no research data or archives needed.
# Usage: Rscript tests/check_ghs_priors.R SOURCE_DIR LIBRARY_DIR OUTPUT_DIR
#    or: Rscript tests/check_ghs_priors.R RUN_ROOT  (source/ and library/ beneath it)
# Optional: CORTREE_REFERENCE_LIBRARY enables exact replay against an older build.
# The source directory must be the source used to build LIBRARY_DIR/CorTree.
args <- commandArgs(trailingOnly = TRUE)
if (!nzchar(Sys.getenv("TZ"))) Sys.setenv(TZ = "UTC")
if (!length(args)) {
  cat("Standalone GHS verification skipped; supply source, installed library, and output directories.\n")
  quit(save = "no", status = 0L)
}
if (length(args) == 1L) {
  root <- normalizePath(args[1L], mustWork = TRUE)
  source_root <- file.path(root, "source")
  library_dir <- file.path(root, "library")
  out <- file.path(root, "verification", "ghs_priors")
} else if (length(args) == 3L) {
  source_root <- args[1L]; library_dir <- args[2L]; out <- args[3L]
} else stop("Usage: Rscript tests/check_ghs_priors.R SOURCE_DIR LIBRARY_DIR OUTPUT_DIR (or RUN_ROOT)")
source_root <- normalizePath(source_root, mustWork = TRUE)
library_dir <- normalizePath(library_dir, mustWork = TRUE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
# Only our previous success marker is removed; a failed rerun cannot look complete.
unlink(file.path(out, "COMPLETE"))
.libPaths(c(library_dir, .libPaths()))
for (dependency in c("Rcpp", "RcppArmadillo", "ape", "mclust")) {
  if (!requireNamespace(dependency, quietly = TRUE)) stop("Required verification dependency is missing: ", dependency)
}
library(CorTree)
stopifnot(normalizePath(find.package("CorTree")) == normalizePath(file.path(library_dir, "CorTree")))
ns <- asNamespace("CorTree")
checks <- list()
record <- function(name, value = TRUE) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(value))
  write.csv(do.call(rbind, checks), file.path(out, "checks.csv"), row.names = FALSE)
  if (!isTRUE(value)) stop("FAILED: ", name)
  cat("PASS", name, "\n")
}
for (sampler in c("CorTree_sampler", "PhyloTree_sampler")) {
  signature <- formals(get(sampler, envir = ns))
  record(paste(sampler, "infinite upper bound default"), identical(signature$ghs_diag_upper, Inf))
  record(paste(sampler, "initial labels remain required"), identical(signature$init_Z, quote(expr = )))
}
# CDF/KS distances below are descriptive calibration tolerances. No IID p-values
# or convergence certification are claimed for correlated Gibbs draws.
cdf_distance <- function(x, cdf, ...) {
  x <- sort(x)
  probability <- cdf(x, ...)
  n <- length(x)
  max(seq_len(n)/n - probability, probability - (seq_len(n) - 1L)/n)
}
empirical_distance <- function(x, y) {
  points <- sort(unique(c(x, y)))
  max(abs(ecdf(x)(points) - ecdf(y)(points)))
}
source_dir <- normalizePath(file.path(source_root, "src"), mustWork = TRUE)
cpp <- paste0(
  '// [[Rcpp::depends(RcppArmadillo)]]\n#include <RcppArmadillo.h>\n',
  '#include "', source_dir, '/PolyaGamma.cpp"\n',
  '#include "', source_dir, '/CorTree.cpp"\n',
  '#include "', source_dir, '/PhyloTree.cpp"\n',
  'template <typename Model> Rcpp::List soft_impl(arma::mat S, int n, double a, double rho, int draws) {\n',
  ' Model model; model.ghs_det_df=a; model.ghs_diag_rate=rho; int p=S.n_rows;\n',
  ' model.ghs_list.emplace_back(p); auto& state=model.ghs_list.front();\n',
  ' arma::mat values(draws,3), inverse_variance(draws,p); double error=0.0, mineig=arma::datum::inf;\n',
  ' for(int i=0;i<draws;++i) { model.GHS_oneSample(S,n,state);\n',
  '  values(i,0)=state.Omega(0,0); values(i,1)=state.Omega(p-1,p-1);\n',
  '  values(i,2)=p>1 ? state.Omega(0,p-1) : 0.0;\n',
  '  inverse_variance.row(i)=(1.0/state.Sigma.diag()).t();\n',
  '  error=std::max(error,arma::abs(state.Omega*state.Sigma-arma::eye(p,p)).max());\n',
  '  mineig=std::min(mineig,arma::eig_sym(state.Omega).min()); }\n',
  ' return Rcpp::List::create(Rcpp::Named("values")=values,\n',
  '  Rcpp::Named("inverse_variance")=inverse_variance,Rcpp::Named("inverse_error")=error,\n',
  '  Rcpp::Named("min_eigenvalue")=mineig); }\n',
  '// [[Rcpp::export]]\n',
  'Rcpp::List actual_soft_draws(arma::mat S,int n,double a,double rho,int draws,bool phylo=false) {\n',
  ' if(phylo) return soft_impl<PhyloTree>(S,n,a,rho,draws);\n',
  ' return soft_impl<CorTree>(S,n,a,rho,draws); }\n',
  '// [[Rcpp::export]]\n',
  'arma::vec actual_component_scales(arma::mat S,arma::mat Q,int n,double shape,double b,int draws) {\n',
  ' arma::vec ans(draws);for(int i=0;i<draws;++i)\n',
  ' ans(i)=cortree::ghs_component_scale(S,Q,n,shape,b);return ans;}\n',
  '// [[Rcpp::export]]\n',
  'arma::vec actual_common_rates(arma::vec t,double shape,double r,double s,int draws) {\n',
  ' arma::vec ans(draws);for(int i=0;i<draws;++i)\n',
  ' ans(i)=cortree::ghs_common_scale_rate(t,shape,r,s);return ans;}\n',
  'template <typename Model> arma::mat hierarchy_prior_impl(int draws) {\n',
  ' Model model;model.ghs_det_df=4;model.ghs_diag_rate=2;\n',
  ' const int K=3; for(int k=0;k<K;++k)model.ghs_list.emplace_back(1);\n',
  ' arma::vec t(K,arma::fill::ones);arma::mat S(1,1,arma::fill::zeros),ans(draws,5);double b=2;\n',
  ' for(int i=0;i<draws;++i) {for(int k=0;k<K;++k){\n',
  '  model.GHS_oneSample(t(k)*S,0,model.ghs_list[k]);\n',
  '  t(k)=cortree::ghs_component_scale(S,model.ghs_list[k].Omega,0,3,b);}\n',
  ' b=cortree::ghs_common_scale_rate(t,3,5,2);\n',
  ' ans(i,0)=model.ghs_list[0].Omega(0,0);ans(i,1)=b;\n',
  ' ans(i,2)=t(0);ans(i,3)=t(1);ans(i,4)=t(0)*model.ghs_list[0].Omega(0,0);}\n',
  ' return ans;}\n',
  '// [[Rcpp::export]]\n',
  'arma::mat actual_hierarchy_prior(int draws,bool phylo=false) {\n',
  ' if(phylo)return hierarchy_prior_impl<PhyloTree>(draws);\n',
  ' return hierarchy_prior_impl<CorTree>(draws);}\n',
  'template <typename Model> Rcpp::List uniform_impl(arma::mat S, int n, double upper, int draws, double det_df) {\n',
  ' Model model; model.ghs_diag_upper=upper; model.ghs_det_df=det_df; int p=S.n_rows;\n',
  ' model.ghs_list.emplace_back(p); auto& state=model.ghs_list.front();\n',
  ' double initial=std::min(1.0,upper/2.0); state.Omega*=initial; state.Sigma/=initial;\n',
  ' arma::mat values(draws,3); double error=0.0, mineig=arma::datum::inf;\n',
  ' for(int i=0;i<draws;++i) { model.GHS_oneSample(S,n,state);\n',
  '  values(i,0)=state.Omega(0,0); values(i,1)=state.Omega(p-1,p-1);\n',
  '  values(i,2)=p>1 ? state.Omega(0,p-1) : 0.0;\n',
  '  error=std::max(error,arma::abs(state.Omega*state.Sigma-arma::eye(p,p)).max());\n',
  '  mineig=std::min(mineig,arma::eig_sym(state.Omega).min()); }\n',
  ' return Rcpp::List::create(Rcpp::Named("values")=values,Rcpp::Named("inverse_error")=error,\n',
  '  Rcpp::Named("min_eigenvalue")=mineig,Rcpp::Named("blocks")=model.uniform_block_updates,\n',
  '  Rcpp::Named("evaluations")=model.uniform_ess_evaluations); }\n',
  '// [[Rcpp::export]]\n',
  'Rcpp::List actual_uniform_draws(arma::mat S,int n,double upper,int draws,bool phylo=false,double det_df=4.0) {\n',
  ' if(phylo) return uniform_impl<PhyloTree>(S,n,upper,draws,det_df);\n',
  ' return uniform_impl<CorTree>(S,n,upper,draws,det_df); }\n',
  '// [[Rcpp::export]]\n',
  'arma::mat noncentral_uniform_block(int draws) {\n',
  ' arma::vec beta(2,arma::fill::zeros), invvar(2), cross(2);\n',
  ' invvar(0)=1.0;invvar(1)=0.7;cross(0)=0.7;cross(1)=-0.4;\n',
  ' arma::mat invA(2,2), result(draws,3);invA(0,0)=1.0;invA(1,1)=1.3;invA(0,1)=invA(1,0)=0.2;\n',
  ' for(int i=0;i<draws;++i) { auto z=cortree::bounded_ghs_block(beta,invA,invvar,cross,0.4,5.5,3.0);\n',
  ' beta=z.beta;result(i,0)=beta(0);result(i,1)=beta(1);result(i,2)=z.gamma; }return result;}\n')
cpp_file <- file.path(out, "actual_ghs_prior_checks.cpp")
writeLines(cpp, cpp_file)
Rcpp::sourceCpp(cpp_file, cacheDir = file.path(out, "cpp_cache"), rebuild = TRUE)

# p=1 reduces exactly to a Gamma law, including empty slots and zero scatter.
scalar_rows <- list()
for (phylo in c(FALSE, TRUE)) {
  name <- if (phylo) "PhyloTree" else "CorTree"
  for (case in list(c(n = 0, S = 0, a = 4, rho = 2),
                   c(n = 7, S = 4, a = 10, rho = 0.25),
                   c(n = 3, S = 0, a = 2, rho = 0.1),
                   c(n = 0, S = 0, a = 0, rho = log(2)),
                   c(n = 0, S = 0, a = 3640.8, rho = 163.836),
                   c(n = 910, S = 81.9, a = 3640.8, rho = 18.204),
                   c(n = 0, S = 0, a = 6518.4, rho = 3259.2))) {
    set.seed(74019 + as.integer(phylo))
    ans <- actual_soft_draws(matrix(case["S"], 1L), as.integer(case["n"]),
      case["a"], case["rho"], 20000L, phylo)
    x <- ans$values[, 1L]
    shape <- (case["n"] + case["a"]) / 2 + 1
    rate <- case["S"] / 2 + case["rho"]
    ks <- cdf_distance(x, pgamma, shape = shape, rate = rate)
    record(paste(name, "scalar Gamma", paste(case, collapse = ":")),
      ks < 0.025 && abs(mean(x) - shape/rate) < 6 * sqrt(shape/rate^2/length(x)) &&
      abs(var(x)/(shape/rate^2) - 1) < 0.12 && ans$inverse_error < 1e-12 &&
      ans$min_eigenvalue > 0)
    scalar_rows[[length(scalar_rows) + 1L]] <- data.frame(sampler = name,
      n = case["n"], scatter = case["S"], a = case["a"], rho = case["rho"],
      mean = mean(x), target_mean = shape/rate, variance = var(x),
      target_variance = shape/rate^2, ks = ks)
  }
}
write.csv(do.call(rbind, scalar_rows), file.path(out, "scalar_checks.csv"), row.names = FALSE)

# Independent p=2 rejection reference for the full joint prior. Diagonal
# Gamma proposals absorb the determinant's diagonal powers, and SPD plus
# (1-off^2/(d1*d2))^(a/2) supply the remaining acceptance factor. This uses
# no production block conditional and also checks the strong prior regime.
prior_rows <- list()
for (prior in list(c(a = 4, rho = 2), c(a = 3640.8, rho = 163.836))) {
  a <- unname(prior["a"]); rho <- unname(prior["rho"])
  set.seed(78517)
  reference <- matrix(numeric(), ncol = 3L)
  batches <- 0L
  while (nrow(reference) < 60000L) {
    batches <- batches + 1L
    if (batches > 2000L) stop("Independent reference sampler exceeded its proposal budget")
    d1 <- rgamma(20000L, a/2 + 1, rate = rho)
    d2 <- rgamma(20000L, a/2 + 1, rate = rho)
    off <- rnorm(20000L) * abs(rcauchy(20000L)) * abs(rcauchy(20000L))
    ratio <- off^2/(d1*d2)
    keep <- is.finite(ratio) & ratio < 1 & runif(20000L) < pmax(0, 1-ratio)^(a/2)
    reference <- rbind(reference, cbind(d1, d2, off)[keep, , drop = FALSE])
  }
  reference <- reference[seq_len(60000L), , drop = FALSE]
  for (phylo in c(FALSE, TRUE)) {
    name <- if (phylo) "PhyloTree" else "CorTree"
    set.seed(38191 + as.integer(phylo))
    ans <- actual_soft_draws(matrix(0, 2L, 2L), 0L, a, rho, 70000L, phylo)
    kept <- ans$values[-seq_len(10000L), , drop = FALSE]
    summaries <- cbind(kept[, 1:2], abs(kept[, 3]), kept[, 3]/sqrt(kept[, 1]*kept[, 2]))
    targets <- cbind(reference[, 1:2], abs(reference[, 3]), reference[, 3]/sqrt(reference[, 1]*reference[, 2]))
    ks <- vapply(1:4, function(j) empirical_distance(summaries[, j], targets[, j]), numeric(1L))
    record(paste(name, "full p2 prior versus independent rejection a=", a),
      max(ks) < 0.05 && ans$inverse_error < 1e-6 && ans$min_eigenvalue > 0)
    gamma <- ans$inverse_variance[-seq_len(10000L), , drop = FALSE]
    gamma_ks <- apply(gamma, 2L, function(x) cdf_distance(x, pgamma, shape = a/2 + 1, rate = rho))
    record(paste(name, "p2 inverse marginal variance Gamma calibration a=", a), max(gamma_ks) < 0.035)
    prior_rows[[length(prior_rows) + 1L]] <- data.frame(sampler = name, a = a, rho = rho,
      statistic = c("diag1", "diag2", "abs_offdiag", "signed_precision_correlation"), ks = ks)
  }
}
write.csv(do.call(rbind, prior_rows), file.path(out, "prior_reference_checks.csv"), row.names = FALSE)

# Exact scalar targets include empty components and positive n with zero scatter.
scalar_rows <- list()
for (phylo in c(FALSE, TRUE)) {
  name <- if (phylo) "PhyloTree" else "CorTree"
  for (case in list(c(n = 0, S = 0, M = 3, a = 0), c(n = 0, S = 0, M = 3, a = 4), c(n = 7, S = 4, M = 3, a = 227.55),
                   c(n = 0, S = 0, M = 100, a = 3640.8), c(n = 200, S = 4, M = 1e-5, a = 6518.4))) {
    set.seed(41951 + as.integer(phylo))
    ans <- actual_uniform_draws(matrix(case["S"], 1L), as.integer(case["n"]), case["M"], 20000L, phylo, case["a"])
    x <- ans$values[, 1L]; shape <- (case["n"] + case["a"]) / 2 + 1; rate <- case["S"] / 2
    cdf <- if (rate == 0) function(q) pmin(1, pmax(0, q / case["M"])^shape) else
      function(q) exp(pgamma(q, shape, rate = rate, log.p = TRUE) - pgamma(case["M"], shape, rate = rate, log.p = TRUE))
    ks <- cdf_distance(x, cdf)
    record(paste(name, "bounded scalar", paste(case, collapse = ":")),
      ks < 0.025 && all(x > 0 & x < case["M"]) && ans$inverse_error < 1e-12)
    scalar_rows[[length(scalar_rows) + 1L]] <- data.frame(sampler = name, n = case["n"], scatter = case["S"], upper = case["M"], ks = ks)
  }
}
write.csv(do.call(rbind, scalar_rows), file.path(out, "bounded_scalar_checks.csv"), row.names = FALSE)

# Independent p=2 reference for a=4: power-law diagonal proposals,
# unchanged off-diagonal hierarchy, acceptance (1 - b^2/(a*d))^2 on SPD.
# This reference does not use a Gibbs sampler or any bounded-block formula.
set.seed(715003)
reference <- matrix(numeric(), ncol = 3L)
batches <- 0L
while (nrow(reference) < 60000L) {
  batches <- batches + 1L
  if (batches > 2000L) stop("Bounded reference sampler exceeded its proposal budget")
  a <- 3 * runif(20000L)^(1 / 3); d <- 3 * runif(20000L)^(1 / 3)
  b <- rnorm(20000L) * abs(rcauchy(20000L)) * abs(rcauchy(20000L))
  keep <- is.finite(b) & b^2 < a * d & runif(20000L) < pmax(0, 1 - b^2 / (a*d))^2
  reference <- rbind(reference, cbind(a, d, b)[keep, , drop = FALSE])
}
reference <- reference[seq_len(60000L), , drop = FALSE]
prior_rows <- list()
for (phylo in c(FALSE, TRUE)) {
  name <- if (phylo) "PhyloTree" else "CorTree"
  set.seed(349613 + as.integer(phylo))
  ans <- actual_uniform_draws(matrix(0, 2L, 2L), 0L, 3, 70000L, phylo)
  kept <- ans$values[-seq_len(10000L), , drop = FALSE]
  summaries <- cbind(kept[, 1:2], abs(kept[, 3]), kept[, 3] / sqrt(kept[, 1] * kept[, 2]))
  targets <- cbind(reference[, 1:2], abs(reference[, 3]), reference[, 3] / sqrt(reference[, 1] * reference[, 2]))
  ks <- vapply(seq_len(ncol(summaries)), function(j) empirical_distance(summaries[, j], targets[, j]), numeric(1L))
  record(paste(name, "p2 full prior versus independent rejection"),
    max(ks) < 0.05 && ans$inverse_error < 1e-6 && ans$min_eigenvalue > 0 &&
    all(kept[, 1:2] > 0 & kept[, 1:2] < 3) && ans$blocks == 140000 && ans$evaluations >= ans$blocks)
  prior_rows[[length(prior_rows) + 1L]] <- data.frame(sampler = name, statistic = c("diag1", "diag2", "abs_offdiag", "signed_correlation"), ks = ks)
}
write.csv(do.call(rbind, prior_rows), file.path(out, "bounded_prior_reference_checks.csv"), row.names = FALSE)

# Nonzero mean and non-diagonal Gaussian precision check ESS centering and Cholesky orientation.
set.seed(34096)
inv_A <- matrix(c(1, 0.2, 0.2, 1.3), 2L)
Q <- 0.4 * inv_A + diag(c(1, 0.7))
C <- solve(Q)
mu <- -C %*% c(0.7, -0.4)
ref_block <- matrix(numeric(), ncol = 3L)
batches <- 0L
while (nrow(ref_block) < 40000L) {
  batches <- batches + 1L
  if (batches > 2000L) stop("Block reference sampler exceeded its proposal budget")
  beta <- sweep(matrix(rnorm(40000L), 20000L, 2L) %*% chol(C), 2L, as.vector(mu), "+")
  gamma <- qgamma(runif(20000L) * pgamma(3, 5.5, rate = 0.2), 5.5, rate = 0.2)
  q <- rowSums(beta * (beta %*% inv_A))
  good <- q + gamma < 3
  ref_block <- rbind(ref_block, cbind(beta, gamma)[good, , drop = FALSE])
}
ref_block <- ref_block[seq_len(40000L), , drop = FALSE]
set.seed(44187)
block <- noncentral_uniform_block(45000L)[-seq_len(5000L), , drop = FALSE]
block_ks <- vapply(1:3, function(j) empirical_distance(block[, j], ref_block[, j]), numeric(1L))
q <- rowSums(block[, 1:2] * (block[, 1:2] %*% inv_A))
record("Noncentral multivariate collapsed block versus independent rejection",
  max(block_ks) < 0.04 && all(q + block[, 3] < 3) &&
  max(abs(cov(block) - cov(ref_block))) < 0.04)
write.csv(data.frame(statistic = c("beta1", "beta2", "gamma"), ks = block_ks), file.path(out, "noncentral_block_checks.csv"), row.names = FALSE)

# Exact scale conditionals exercise the shared production helpers. A p=2
# scatter with cross terms detects using only diagonals in tr(SQ).
Q <- matrix(c(2, 0.4, 0.4, 1.5), 2L)
S <- matrix(c(4, -0.7, -0.7, 3), 2L)
set.seed(31191)
scales <- actual_component_scales(S, Q, 7L, 3, 2, 20000L)
shape <- 3 + 7 * 2 / 2
rate <- 2 + sum(S * Q) / 2
record("Component scale conditional uses n*p/2 and full tr(SQ)",
  cdf_distance(scales, pgamma, shape = shape, rate = rate) < 0.025)
set.seed(31192)
empty_scales <- actual_component_scales(matrix(0, 2L, 2L), Q, 0L, 3, 2, 20000L)
record("Empty component scale conditional remains proper",
  cdf_distance(empty_scales, pgamma, shape = 3, rate = 2) < 0.025)
set.seed(31193)
rates <- actual_common_rates(c(0.2, 1, 2, 4), 3, 2, 1, 20000L)
record("Shared rate conditional includes all K component scales",
  cdf_distance(rates, pgamma, shape = 2 + 4*3, rate = 1 + 7.2) < 0.025)

# Under the joint prior b~Gamma(r,s), t/s~BetaPrime(c,r), and independent
# components have correlation c/(c+r-1). Check the coupled chain, not just
# isolated Gamma draws. r=5 gives finite moments suitable for this test.
hierarchy_rows <- list()
for (phylo in c(FALSE, TRUE)) {
  name <- if (phylo) "PhyloTree" else "CorTree"
  set.seed(21367 + as.integer(phylo))
  kept <- actual_hierarchy_prior(70000L, phylo)[-seq_len(10000L), , drop = FALSE]
  beta_prime_cdf <- function(x) pbeta(x/(2+x), 3, 5)
  ks <- c(Q = cdf_distance(kept[, 1L], pgamma, shape = 3, rate = 2),
    b = cdf_distance(kept[, 2L], pgamma, shape = 5, rate = 2),
    t = cdf_distance(kept[, 3L], beta_prime_cdf))
  shared_cor <- cor(kept[, 3L], kept[, 4L])
  record(paste(name, "joint hierarchy prior and induced scale dependence"),
    max(ks) < 0.035 && abs(shared_cor - 3/7) < 0.07 &&
    abs(cor(kept[, 1L], kept[, 3L])) < 0.04 &&
    all(is.finite(kept)) && all(kept > 0))
  hierarchy_rows[[length(hierarchy_rows) + 1L]] <- data.frame(sampler = name,
    ks_Q = ks[1L], ks_b = ks[2L], ks_t = ks[3L], component_scale_correlation = shared_cor,
    target_scale_correlation = 3/7, mean_effective_precision = mean(kept[, 5L]),
    target_mean_effective_precision = 1.5*1.5)
}
write.csv(do.call(rbind, hierarchy_rows), file.path(out, "hierarchy_prior_checks.csv"), row.names = FALSE)

# Small integrated fits exercise effective precision tQ, empty-slot refresh,
# retained scale traces, proper-path safeguards, and both tree samplers.
set.seed(126)
X <- matrix(rpois(4L*8L, 4), 4L, 8L)
tree <- ape::stree(8L, type = "balanced"); colnames(X) <- tree$tip.label
cor_args <- list(X = X, n_clus = 6L, tree_depth = 3L, cutoff_layer = 1L,
  total_iter = 16L, burnin = 4L, warm_start = 0L, init_Z = rep(0L, 4L),
  c_sigma2_vec = 10, sigma_mu2 = 0.1, cov_interval = 2L,
  ghs_diag_rate = 2, ghs_det_df = 4)
phy_args <- cor_args; phy_args$X <- phy_args$tree_depth <- NULL
phy_args$count_data <- X; phy_args$tree <- tree; phy_args$save_sigma_inv_trace <- TRUE
for (sampler in c("CorTree_sampler", "PhyloTree_sampler")) {
  fun <- get(sampler, envir = ns)
  pars <- if (sampler == "CorTree_sampler") cor_args else phy_args
  for (hierarchical in c(FALSE, TRUE)) {
    pars$ghs_scale_hierarchy <- hierarchical
    set.seed(52918)
    m <- do.call(fun, pars)$mcmc
    record(paste(sampler, "proper metadata and no safeguards", hierarchical),
      isTRUE(m$ghs_prior$proper) && isTRUE(m$ghs_prior$active) &&
      m$ghs_prior$det_df == 4 && m$ghs_prior$diag_rate == 2 &&
      !m$covariance_safeguards$enabled && sum(m$covariance_safeguards$regularizations) == 0 &&
      sum(m$covariance_safeguards$updates_skipped) == 0 && all(is.finite(m$loglik)))
    record(paste(sampler, "empty precision refresh", hierarchical), sum(m$empty_precision_updates) > 0)
    for (cube in m$Sigma_inv) for (k in seq_len(dim(cube)[3L])) chol(cube[, , k])
    record(paste(sampler, "all saved effective precisions SPD", hierarchical))
    if (hierarchical) {
      h <- m$ghs_scale
      record(paste(sampler, "scale trace shape and values"),
        isTRUE(h$active) && identical(dim(h$t), c(6L, 12L)) && length(h$b) == 12L &&
        all(is.finite(h$t)) && all(h$t > 0) && all(is.finite(h$b)) && all(h$b > 0) &&
        h$shape == 3 && h$rate_shape == 2 && h$rate_rate == 1)
      record(paste(sampler, "scale updates follow covariance interval"),
        identical(h$t[, seq(1L, 11L, by = 2L)], h$t[, seq(2L, 12L, by = 2L)]) &&
        identical(as.vector(h$b)[seq(1L, 11L, by = 2L)], as.vector(h$b)[seq(2L, 12L, by = 2L)]))
      # K>N guarantees two slots are empty at every update. Their t and Q must
      # still refresh; at least one never-occupied slot makes this observable.
      unused <- setdiff(0:5, as.integer(m$Z)) + 1L
      record(paste(sampler, "never occupied slot scale refresh"),
        length(unused) > 0 && all(apply(h$t[unused, , drop = FALSE], 1L, function(x) length(unique(x))) > 1L))
    }
  }
  for (bad in list(list(ghs_diag_rate = 0), list(ghs_diag_upper = 3),
                   list(ghs_scale_shape = 0), list(ghs_scale_shape = NA_real_),
                   list(ghs_scale_rate_shape = -1), list(ghs_scale_rate_rate = Inf))) {
    invalid <- modifyList(pars, bad)
    error <- tryCatch({do.call(fun, invalid); ""}, error = conditionMessage)
    record(paste(sampler, "invalid hierarchy rejected", paste(names(bad), collapse = ":")),
      nzchar(error) && grepl("ghs_|GHS", error))
  }
  independent <- pars; independent$all_ind <- TRUE; independent$ghs_scale_hierarchy <- FALSE
  set.seed(8112); without <- do.call(fun, independent)$mcmc
  independent$ghs_scale_hierarchy <- TRUE
  set.seed(8112); with <- do.call(fun, independent)$mcmc
  record(paste(sampler, "hierarchy inactive for independent mode"),
    identical(without$Z, with$Z) && identical(without$loglik, with$loglik) &&
    !with$ghs_scale$active && is.null(with$ghs_scale$t) && is.null(with$ghs_scale$b))
}

# Bounded modes must update empty slots, remain inside open diagonal support,
# and avoid the legacy flat-prior numerical safeguards.
for (sampler in c("CorTree_sampler", "PhyloTree_sampler")) {
  fun <- get(sampler, envir = ns)
  pars <- if (sampler == "CorTree_sampler") cor_args else phy_args
  for (a in c(0, 4, 3640.8)) {
    bounded_args <- modifyList(pars, list(ghs_diag_rate = 0, ghs_diag_upper = 0.5,
      ghs_det_df = a, ghs_scale_hierarchy = FALSE))
    set.seed(592)
    m <- do.call(fun, bounded_args)$mcmc
    record(paste(sampler, "bounded empty updates and proper metadata a=", a),
      sum(m$empty_precision_updates) > 0 && isTRUE(m$ghs_prior$proper) &&
      isTRUE(m$ghs_prior$active) && m$ghs_prior$diag_upper == 0.5 &&
      m$ghs_prior$det_df == a && !m$covariance_safeguards$enabled &&
      sum(m$covariance_safeguards$regularizations) == 0 &&
      sum(m$covariance_safeguards$updates_skipped) == 0 && all(is.finite(m$loglik)) &&
      m$ghs_uniform$block_updates > 0 && m$ghs_uniform$ess_evaluations >= m$ghs_uniform$block_updates)
    for (cube in m$Sigma_inv) for (k in seq_len(dim(cube)[3L])) {
      chol(cube[, , k])
      stopifnot(all(is.finite(cube[, , k])), all(diag(cube[, , k]) > 0 & diag(cube[, , k]) < 0.5))
    }
    record(paste(sampler, "all retained bounded precisions in support a=", a))
  }
  invalid_options <- list(list(ghs_diag_rate = -1), list(ghs_diag_rate = Inf),
    list(ghs_diag_upper = 0), list(ghs_diag_upper = NA_real_),
    list(ghs_diag_rate = 1, ghs_diag_upper = 3),
    list(ghs_diag_rate = 0, ghs_diag_upper = Inf, ghs_det_df = 1),
    list(ghs_det_df = -1), list(ghs_det_df = Inf),
    list(ghs_diag_rate = 0, ghs_diag_upper = 3, warm_start = 2L))
  for (bad in invalid_options) {
    error <- tryCatch({do.call(fun, modifyList(pars, bad)); ""}, error = conditionMessage)
    record(paste(sampler, "invalid prior options", paste(names(bad), unlist(bad), collapse = ":")),
      nzchar(error) && grepl("ghs_|GHS", error))
  }
  independent <- modifyList(pars, list(all_ind = TRUE, ghs_diag_rate = 0,
    ghs_diag_upper = Inf, ghs_det_df = 0, ghs_scale_hierarchy = FALSE))
  set.seed(8112); baseline <- do.call(fun, independent)$mcmc
  modes <- list(exponential = list(ghs_diag_rate = log(2)),
    uniform = list(ghs_diag_upper = 3),
    bounded_determinant = list(ghs_diag_upper = 3, ghs_det_df = 3640.8),
    soft_determinant = list(ghs_diag_rate = 163.836, ghs_det_df = 3640.8),
    hierarchy = list(ghs_diag_rate = 2, ghs_det_df = 4, ghs_scale_hierarchy = TRUE))
  for (mode in names(modes)) {
    set.seed(8112); observed <- do.call(fun, modifyList(independent, modes[[mode]]))$mcmc
    # Diagnostic metadata legitimately changes; posterior traces must not.
    keep <- intersect(c("mu", "phi", "sigma2_vec", "pi", "Z", "loglik"), names(baseline))
    record(paste(sampler, "all_ind exact trace parity", mode),
      identical(baseline[keep], observed[keep]) && !observed$ghs_prior$active &&
      sum(observed$empty_precision_updates) == 0)
  }
}

# Wrapper checks include user-facing multistart and held-out refit entry points.
multi_args <- cor_args; multi_args$n_clus <- 2L; multi_args$n_start <- 1L
multi_args$init_Z <- NULL; multi_args$init_Z_list <- list(c(0L, 0L, 1L, 1L))
multi_args$discard_collapsed_at_burnin <- FALSE; multi_args$verbose <- FALSE
multi_args$ghs_scale_hierarchy <- TRUE
multi_args$ghs_scale_shape <- 4; multi_args$ghs_scale_rate_shape <- 5; multi_args$ghs_scale_rate_rate <- 2
set.seed(612)
multi <- do.call(get("CorTree_sampler_randominit", envir = ns), multi_args)
record("Random initialization forwards hierarchy parameters", multi$all_fit[[1L]]$mcmc$ghs_scale$shape == 4 &&
  multi$all_fit[[1L]]$mcmc$ghs_scale$rate_shape == 5 && multi$all_fit[[1L]]$mcmc$ghs_scale$rate_rate == 2)
held_args <- multi_args; held_args$train_idx <- 1:2; held_args$test_idx <- 3:4
held_args$refit_full_data <- TRUE
set.seed(613)
held <- do.call(get("CorTree_sampler_randominit_heldout", envir = ns), held_args)
record("Held-out training and refit forward hierarchy parameters",
  held$best_fit_heldout$mcmc$ghs_scale$shape == 4 &&
  held$full_fit_from_best_start_init$mcmc$ghs_scale$shape == 4)

# The other proper modes must also reach both multistart wrapper entry points.
for (mode in list(list(ghs_diag_rate = log(2), ghs_diag_upper = Inf, ghs_det_df = 0),
                  list(ghs_diag_rate = 0, ghs_diag_upper = 3, ghs_det_df = 0),
                  list(ghs_diag_rate = 0, ghs_diag_upper = 3, ghs_det_df = 4),
                  list(ghs_diag_rate = 163.836, ghs_diag_upper = Inf, ghs_det_df = 3640.8))) {
  p <- modifyList(multi_args, c(mode, list(ghs_scale_hierarchy = FALSE)))
  set.seed(614); fit <- do.call(get("CorTree_sampler_randominit", envir = ns), p)
  got <- fit$all_fit[[1L]]$mcmc$ghs_prior
  record(paste("Multistart forwards proper prior", paste(unlist(mode), collapse = ":")),
    got$diag_rate == mode$ghs_diag_rate && got$diag_upper == mode$ghs_diag_upper && got$det_df == mode$ghs_det_df)
  p$train_idx <- 1:2; p$test_idx <- 3:4; p$refit_full_data <- TRUE
  set.seed(615); fit <- do.call(get("CorTree_sampler_randominit_heldout", envir = ns), p)
  traces <- list(fit$best_fit_heldout$mcmc$ghs_prior, fit$full_fit_from_best_start_init$mcmc$ghs_prior)
  record(paste("Held-out and refit forward proper prior", paste(unlist(mode), collapse = ":")),
    all(vapply(traces, function(got) got$diag_rate == mode$ghs_diag_rate &&
      got$diag_upper == mode$ghs_diag_upper && got$det_df == mode$ghs_det_df, logical(1L))))
}

# Existing priors must reproduce their previous random draws exactly. Separate
# R processes prevent accidental linkage to the wrong package's native code.
baseline_lib <- Sys.getenv("CORTREE_REFERENCE_LIBRARY", "")
if (nzchar(baseline_lib)) {
baseline_lib <- normalizePath(baseline_lib, mustWork = TRUE)
stopifnot(dir.exists(file.path(baseline_lib, "CorTree")))
parity_script <- file.path(out, "prior_parity.R")
writeLines(c(
  'a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); library(CorTree)',
  'set.seed(491); X <- matrix(rpois(8*8,3),8,8); tree <- ape::stree(8,type="balanced"); colnames(X) <- tree$tip.label',
  'p <- list(n_clus=3L,cutoff_layer=1L,total_iter=12L,burnin=4L,init_Z=rep(0:2,length.out=8),cov_interval=2L)',
  'result <- list(); for(case in 1:4) { p$ghs_diag_rate <- c(0,log(2),0,0)[case]; p$ghs_diag_upper <- if(case>=3) 3 else Inf; p$ghs_det_df <- if(case==4) 4 else 0',
  ' set.seed(517); co <- do.call(CorTree_sampler,c(list(X=X,tree_depth=3L),p))$mcmc',
  ' set.seed(617); ph <- do.call(PhyloTree_sampler,c(list(count_data=X,tree=tree,save_sigma_inv_trace=TRUE),p))$mcmc',
  ' result[[length(result)+1]] <- list(co=co,ph=ph) }; saveRDS(result,a[2])'
), parity_script)
for (name in c("baseline", "updated")) {
  lib <- if (name == "baseline") baseline_lib else library_dir
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(parity_script), shQuote(lib), shQuote(file.path(out, paste0(name, "_parity.rds")))),
    stdout = file.path(out, paste0(name, "_parity.log")), stderr = file.path(out, paste0(name, "_parity.log")))
  stopifnot(status == 0L)
}
old <- readRDS(file.path(out, "baseline_parity.rds"))
new <- readRDS(file.path(out, "updated_parity.rds"))
for (i in seq_along(old)) for (name in c("co", "ph")) {
  target <- old[[i]][[name]]; observed <- new[[i]][[name]][names(target)]
  observed$ghs_prior <- observed$ghs_prior[names(target$ghs_prior)]
  record(paste(name, "exact prior replay", c("flat", "exponential", "uniform", "bounded determinant")[i]),
    identical(target, observed))
}
} else {
  cat("SKIP optional archived-build replay: CORTREE_REFERENCE_LIBRARY is not set.\n")
}
write.csv(do.call(rbind, checks), file.path(out, "checks.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(out, "COMPLETE"))
cat("All standalone GHS prior verification checks passed.\n")
