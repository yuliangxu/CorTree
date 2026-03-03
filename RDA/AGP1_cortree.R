rm(list = ls())
devtools::load_all()
set.seed(2026)

project_dir <- normalizePath(".")
data_path <- file.path(project_dir, "data", "Cluster_AG_subsample.rds")
out_dir <- file.path(project_dir, "data")
out_rds <- file.path(out_dir, "AGP1_phylotree_fit.rds")


ag_obj <- readRDS(data_path)
count_data <- t(as.matrix(ag_obj$otu_top))



ntile_base <- function(x, n) {
  as.integer(cut(rank(x, ties.method = "first"), breaks = n, labels = FALSE))
}

n_clus <- 5L
init_Z <- ntile_base(rowSums(count_data), n_clus) - 1L

tree_info <- aggregate_tree_counts(count_data, ag_obj$tree)
tree_depth <- as.integer(max(tree_info$depth))

cutoff_layer <- min(3L, tree_depth)
burnin <- 5
total_iter <- burnin + 5
warm_start <- 0
c_sigma2_vec <- 3
sigma_mu2 <- 0.5
all_ind <- FALSE
cov_interval <- 3L

fit <- PhyloTree_sampler(
  count_data = count_data,
  tree = ag_obj$tree,
  n_clus = n_clus,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  warm_start = warm_start,
  init_Z = init_Z,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  all_ind = all_ind,
  cov_interval = cov_interval
)

Z_post <- fit$mcmc$Z
cluster_hat <- apply(Z_post, 1, function(z) {
  as.integer(names(which.max(table(z))))
})

summary_obj <- list(
  elapsed_sec = fit$elapsed,
  n_sample = nrow(count_data),
  n_taxa = ncol(count_data),
  n_clus = n_clus,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  cluster_size = table(cluster_hat),
  loglik_tail_mean = mean(tail(fit$mcmc$loglik, n = max(1L, floor(total_iter * 0.2))))
)

result <- list(
  input = list(
    data_path = data_path,
    keep_sample = keep_sample,
    sampler_args = list(
      n_clus = n_clus,
      cutoff_layer = cutoff_layer,
      total_iter = total_iter,
      burnin = burnin,
      warm_start = warm_start,
      c_sigma2_vec = c_sigma2_vec,
      sigma_mu2 = sigma_mu2,
      all_ind = all_ind,
      cov_interval = cov_interval
    )
  ),
  fit = fit,
  summary = summary_obj,
  cluster_hat = cluster_hat
)

saveRDS(result, out_rds)

pdf(out_pdf, width = 9, height = 4)
par(mfrow = c(1, 2))
plot(fit$mcmc$loglik, type = "l", xlab = "Iteration", ylab = "Log-likelihood", main = "PhyloTree log-likelihood")
matplot(t(fit$mcmc$pi), type = "l", lty = 1, xlab = "Saved draw", ylab = "pi", main = "Mixture proportions")
dev.off()

cat("PhyloTree run complete\n")
cat("Elapsed (sec):", fit$elapsed, "\n")
cat("Saved results:", out_rds, "\n")
cat("Saved diagnostics:", out_pdf, "\n")
print(summary_obj$cluster_size)
