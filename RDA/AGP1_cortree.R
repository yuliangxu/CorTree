# Sys.getenv("PATH")
# Sys.getenv("LD_LIBRARY_PATH")
# system("which g++")
# system("g++ --version")

devtools::load_all()


rm(list = ls())
devtools::load_all()
set.seed(2026)

project_dir <- normalizePath(".")
data_path <- file.path(project_dir, "data", "Cluster_AG_subsample.rds")
out_dir = "/cwork/yx306/CorTree"
out_rds <- file.path(out_dir, "AGP1_phylotree_fit.rds")
out_pdf <- file.path(out_dir, "AGP1_phylotree_diagnostics.pdf")
out_cov_assoc <- file.path(out_dir, "AGP1_cluster_covariate_association.csv")


ag_obj <- readRDS(data_path)
count_data <- t(as.matrix(ag_obj$otu_top))
keep_sample <- rownames(count_data)
ag_cov <- ag_obj$ag_fecal

n_clus <- 10L
init_Z <- ntile_base(rowSums(count_data), n_clus) - 1L

tree_info <- aggregate_tree_counts(count_data, ag_obj$tree)
tree_depth <- as.integer(max(tree_info$depth))

cutoff_layer <- min(3L, tree_depth)
burnin <- 150
total_iter <- burnin + 50
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
saveRDS(fit, file.path(out_dir, "AGP1_phylotree_fit.rds"))

fit_size_gb <- as.numeric(object.size(fit)) / 1024^3; print(fit_size_gb)
elapsed_min <- fit$elapsed / 60; print(elapsed_min)
pi_trace_info <- get_pi_trace_plot_data(
  fit = fit,
  sampler_args = list(
    total_iter = total_iter,
    burnin = burnin
  )
)

Z_post <- fit$mcmc$Z
cluster_hat <- apply(Z_post, 1, function(z) {
  as.integer(names(which.max(table(z))))
})

umap_input <- log1p(count_data)
umap_neighbors <- max(5L, min(30L, nrow(umap_input) - 1L))
umap_embedding <- uwot::umap(
  umap_input,
  n_neighbors = umap_neighbors,
  min_dist = 0.2,
  metric = "euclidean",
  scale = TRUE,
  verbose = FALSE,
  ret_model = FALSE
)
colnames(umap_embedding) <- c("UMAP1", "UMAP2")
rownames(umap_embedding) <- rownames(count_data)

covariate_assoc <- analyze_covariate_association(
  cluster = cluster_hat,
  covariates = ag_cov,
  sample_ids = keep_sample
)

top_covariate <- if (nrow(covariate_assoc) > 0L) covariate_assoc[1, , drop = FALSE] else NULL

summary_obj <- list(
  elapsed_sec = fit$elapsed,
  elapsed_min = elapsed_min,
  fit_size_gb = fit_size_gb,
  n_sample = nrow(count_data),
  n_taxa = ncol(count_data),
  n_clus = n_clus,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  cluster_size = table(cluster_hat),
  loglik_tail_mean = mean(tail(fit$mcmc$loglik, n = max(1L, floor(total_iter * 0.2)))),
  strongest_covariate = top_covariate
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
  cluster_hat = cluster_hat,
  umap = umap_embedding,
  covariate_association = covariate_assoc
)

saveRDS(result, out_rds)
utils::write.csv(covariate_assoc, out_cov_assoc, row.names = FALSE)

pdf(out_pdf, width = 10, height = 8)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))
par(mar = c(4, 4, 3, 1))
plot(fit$mcmc$loglik, type = "l", xlab = "Iteration", ylab = "Log-likelihood", main = "PhyloTree log-likelihood")
matplot(
  x = pi_trace_info$iteration_index,
  y = t(pi_trace_info$pi_chain),
  type = "l",
  lty = 1,
  xlab = "MCMC iteration",
  ylab = "pi",
  main = if (identical(pi_trace_info$trace_source, "full")) {
    "Mixture proportions (full trace)"
  } else {
    "Mixture proportions (saved trace)"
  }
)
plot(
  umap_embedding[, 1],
  umap_embedding[, 2],
  col = cluster_hat + 1L,
  pch = 19,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = "UMAP of log1p(count_data)"
)
legend(
  "topright",
  legend = paste("Cluster", sort(unique(cluster_hat))),
  col = sort(unique(cluster_hat)) + 1L,
  pch = 19,
  bty = "n",
  cex = 0.8
)
dev.off()

cat("PhyloTree run complete\n")
cat("Elapsed (sec):", fit$elapsed, "\n")
cat("Elapsed (min):", elapsed_min, "\n")
cat("Fit size (GB):", fit_size_gb, "\n")
cat("Saved results:", out_rds, "\n")
cat("Saved diagnostics:", out_pdf, "\n")
cat("Saved covariate associations:", out_cov_assoc, "\n")
print(summary_obj$cluster_size)
if (!is.null(top_covariate)) {
  cat("Strongest covariate association:\n")
  print(top_covariate)
}
