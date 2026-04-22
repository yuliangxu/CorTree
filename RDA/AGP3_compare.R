rm(list = ls())

devtools::load_all()
set.seed(2026)

project_dir <- normalizePath(".")
out_dir <- "/cwork/yx306/CorTree"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(project_dir, "data", "Cluster_AG_subsample.rds")
cortree_bundle_candidates <- unique(c(
  file.path(out_dir, "AGP1_phylotree_fit.rds")
))
cortree_bundle_candidates <- cortree_bundle_candidates[file.exists(cortree_bundle_candidates)]
if (length(cortree_bundle_candidates) == 0L) {
  stop(
    "Could not find AGP1 CorTree result bundle: AGP1_phylotree_fit.rds. ",
    "Please run RDA/AGP1_cortree.R first."
  )
}
cortree_bundle_path <- cortree_bundle_candidates[[1L]]

out_bundle <- file.path(out_dir, "AGP3_compare_bundle.rds")
out_cluster_csv <- file.path(out_dir, "AGP3_compare_cluster_assignments.csv")
out_timing_csv <- file.path(out_dir, "AGP3_compare_timing.csv")
out_method_summary_csv <- file.path(out_dir, "AGP3_compare_method_summary.csv")
out_crosstab_csv <- file.path(out_dir, "AGP3_compare_crosstab_vs_cortree.csv")
out_cov_assoc_csv <- file.path(out_dir, "AGP3_compare_covariate_association.csv")
out_cov_assoc_rds <- file.path(out_dir, "AGP3_compare_covariate_association.rds")
out_pdf <- file.path(out_dir, "AGP3_compare_umap.pdf")
out_dmm_fit <- file.path(out_dir, "AGP3_dmm_fit.rds")
out_indtree_fit <- file.path(out_dir, "AGP3_indtree_fit.rds")
out_dmm_posterior <- file.path(out_dir, "AGP3_dmm_posterior.rds")
out_indtree_pi <- file.path(out_dir, "AGP3_indtree_pi_trace_postburnin.rds")

if (!file.exists(data_path)) {
  stop("Could not find AGP data: ", data_path)
}

if (!requireNamespace("DirichletMultinomial", quietly = TRUE)) {
  stop(
    "AGP3_compare.R needs package 'DirichletMultinomial'. ",
    "Please install it before running the comparison."
  )
}
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop(
    "AGP3_compare.R needs package 'uwot'. ",
    "Please install it before running the comparison."
  )
}

extract_mode_cluster <- function(Z_chain) {
  apply(Z_chain, 1, function(z) {
    as.integer(names(which.max(table(z))))
  })
}

extract_dmm_posterior <- function(dmm_fit) {
  posterior <- tryCatch(
    as.matrix(DirichletMultinomial::mixture(dmm_fit)),
    error = function(e) NULL
  )
  if (is.null(posterior)) {
    stop("Could not extract per-sample posterior cluster probabilities from the DMM fit.")
  }

  if (nrow(posterior) < ncol(posterior)) {
    posterior <- t(posterior)
  }
  storage.mode(posterior) <- "numeric"
  posterior
}

cluster_size_df <- function(cluster, method) {
  tab <- table(cluster)
  data.frame(
    method = method,
    cluster = as.integer(names(tab)),
    size = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

ag_obj <- readRDS(data_path)
count_data <- t(as.matrix(ag_obj$otu_top))
keep_sample <- rownames(count_data)
ag_cov <- ag_obj$ag_fecal

cov_clean_info <- clean_covariate_dataframe(ag_cov)
ag_cov_clean <- cov_clean_info$covariates

cortree_bundle <- readRDS(cortree_bundle_path)
if (is.null(cortree_bundle$fit)) {
  stop("AGP1 bundle does not contain $fit: ", cortree_bundle_path)
}
cortree_fit <- cortree_bundle$fit
cortree_sampler_args <- cortree_bundle$input$sampler_args

if (!is.null(cortree_bundle$cluster_hat)) {
  cluster_cortree <- as.integer(cortree_bundle$cluster_hat)
} else if (!is.null(cortree_fit$mcmc$Z)) {
  cluster_cortree <- extract_mode_cluster(cortree_fit$mcmc$Z)
} else {
  stop("Could not recover CorTree hard cluster assignments from AGP1 bundle.")
}

if (length(cluster_cortree) != nrow(count_data)) {
  stop("CorTree cluster assignment length does not match sample count.")
}

if (!is.null(cortree_bundle$umap)) {
  umap_embedding <- as.matrix(cortree_bundle$umap)
} else {
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
}

n_clus <- if (!is.null(cortree_sampler_args$n_clus)) {
  as.integer(cortree_sampler_args$n_clus)
} else {
  length(unique(cluster_cortree))
}
init_Z <- ntile_base(rowSums(count_data), n_clus) - 1L

tree_info <- aggregate_tree_counts(count_data, ag_obj$tree)
tree_depth <- as.integer(max(tree_info$depth))
cutoff_layer <- if (!is.null(cortree_sampler_args$cutoff_layer)) {
  as.integer(cortree_sampler_args$cutoff_layer)
} else {
  min(3L, tree_depth)
}
burnin <- if (!is.null(cortree_sampler_args$burnin)) {
  as.integer(cortree_sampler_args$burnin)
} else {
  150L
}
total_iter <- if (!is.null(cortree_sampler_args$total_iter)) {
  as.integer(cortree_sampler_args$total_iter)
} else {
  burnin + 50L
}
warm_start <- if (!is.null(cortree_sampler_args$warm_start)) {
  as.integer(cortree_sampler_args$warm_start)
} else {
  0L
}
c_sigma2_vec <- if (!is.null(cortree_sampler_args$c_sigma2_vec)) {
  cortree_sampler_args$c_sigma2_vec
} else {
  3
}
sigma_mu2 <- if (!is.null(cortree_sampler_args$sigma_mu2)) {
  cortree_sampler_args$sigma_mu2
} else {
  0.5
}
cov_interval <- if (!is.null(cortree_sampler_args$cov_interval)) {
  as.integer(cortree_sampler_args$cov_interval)
} else {
  3L
}

# DMM ------------------------------------------------
dmm_time <- system.time({
  dmm_fit <- DirichletMultinomial::dmn(count_data, k = n_clus)
})
dmm_posterior <- extract_dmm_posterior(dmm_fit)
if (nrow(dmm_posterior) != nrow(count_data)) {
  stop("DMM posterior row count does not match the number of samples.")
}
colnames(dmm_posterior) <- paste0("cluster_", seq_len(ncol(dmm_posterior)) - 1L)
rownames(dmm_posterior) <- rownames(count_data)
cluster_dmm <- max.col(dmm_posterior, ties.method = "first") - 1L
dmm_mixture_weight <- tryCatch(
  as.numeric(DirichletMultinomial::mixturewt(dmm_fit)),
  error = function(e) rep(NA_real_, n_clus)
)

# IndTree ------------------------------------------------
indtree_time <- system.time({
  indtree_fit <- PhyloTree_sampler(
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
    cov_interval = cov_interval,
    all_ind = TRUE
  )
})
cluster_indtree <- extract_mode_cluster(indtree_fit$mcmc$Z)
indtree_pi_trace <- extract_pi_trace(indtree_fit, burnin = burnin)

covariate_assoc_list <- list(
  cortree = analyze_covariate_association(cluster_cortree, ag_cov_clean, keep_sample),
  dmm = analyze_covariate_association(cluster_dmm, ag_cov_clean, keep_sample),
  indtree = analyze_covariate_association(cluster_indtree, ag_cov_clean, keep_sample)
)
covariate_assoc_df <- do.call(
  rbind,
  lapply(names(covariate_assoc_list), function(method_name) {
    assoc_df <- covariate_assoc_list[[method_name]]
    if (nrow(assoc_df) == 0L) {
      return(NULL)
    }
    cbind(method = method_name, assoc_df, stringsAsFactors = FALSE)
  })
)
if (is.null(covariate_assoc_df)) {
  covariate_assoc_df <- data.frame()
}

cluster_assignments <- data.frame(
  SampleID = keep_sample,
  CorTree = cluster_cortree,
  DMM = cluster_dmm,
  IndTree = cluster_indtree,
  UMAP1 = umap_embedding[, 1],
  UMAP2 = umap_embedding[, 2],
  stringsAsFactors = FALSE
)

timing_df <- data.frame(
  method = c("CorTree", "DMM", "IndTree"),
  elapsed_sec = c(
    as.numeric(cortree_fit$elapsed),
    unname(dmm_time[["elapsed"]]),
    as.numeric(if (!is.null(indtree_fit$elapsed)) indtree_fit$elapsed else indtree_time[["elapsed"]])
  ),
  elapsed_min = c(
    as.numeric(cortree_fit$elapsed) / 60,
    unname(dmm_time[["elapsed"]]) / 60,
    as.numeric(if (!is.null(indtree_fit$elapsed)) indtree_fit$elapsed else indtree_time[["elapsed"]]) / 60
  ),
  stringsAsFactors = FALSE
)
timing_df$speed_vs_cortree <- timing_df$elapsed_sec / timing_df$elapsed_sec[timing_df$method == "CorTree"]

method_summary <- data.frame(
  method = c("CorTree", "DMM", "IndTree"),
  ari_vs_cortree = c(
    1,
    adjusted_rand_index(cluster_cortree, cluster_dmm),
    adjusted_rand_index(cluster_cortree, cluster_indtree)
  ),
  n_clusters_used = c(
    length(unique(cluster_cortree)),
    length(unique(cluster_dmm)),
    length(unique(cluster_indtree))
  ),
  loglik = c(
    if (!is.null(cortree_fit$mcmc$loglik)) tail(cortree_fit$mcmc$loglik, 1L) else NA_real_,
    tryCatch(as.numeric(stats::logLik(dmm_fit)), error = function(e) NA_real_),
    if (!is.null(indtree_fit$mcmc$loglik)) tail(indtree_fit$mcmc$loglik, 1L) else NA_real_
  ),
  stringsAsFactors = FALSE
)

crosstab_vs_cortree <- rbind(
  transform(
    as.data.frame(table(CorTree = cluster_cortree, Other = cluster_dmm)),
    method = "DMM"
  ),
  transform(
    as.data.frame(table(CorTree = cluster_cortree, Other = cluster_indtree)),
    method = "IndTree"
  )
)
names(crosstab_vs_cortree) <- c("cortree_cluster", "other_cluster", "n", "method")
crosstab_vs_cortree <- crosstab_vs_cortree[, c("method", "cortree_cluster", "other_cluster", "n")]

cluster_sizes <- do.call(
  rbind,
  list(
    cluster_size_df(cluster_cortree, "CorTree"),
    cluster_size_df(cluster_dmm, "DMM"),
    cluster_size_df(cluster_indtree, "IndTree")
  )
)

comparison_bundle <- list(
  input = list(
    data_path = data_path,
    cortree_bundle_path = cortree_bundle_path,
    keep_sample = keep_sample,
    sampler_args = list(
      n_clus = n_clus,
      cutoff_layer = cutoff_layer,
      total_iter = total_iter,
      burnin = burnin,
      warm_start = warm_start,
      c_sigma2_vec = c_sigma2_vec,
      sigma_mu2 = sigma_mu2,
      cov_interval = cov_interval
    )
  ),
  cluster_assignments = cluster_assignments,
  cluster_sizes = cluster_sizes,
  timing = timing_df,
  method_summary = method_summary,
  crosstab_vs_cortree = crosstab_vs_cortree,
  covariate_association = covariate_assoc_list,
  normalization_report = cov_clean_info$normalization_report,
  umap = umap_embedding,
  posterior = list(
    dmm = dmm_posterior
  ),
  trace = list(
    indtree_pi_postburnin = indtree_pi_trace
  ),
  fit = list(
    cortree = cortree_fit,
    dmm = dmm_fit,
    indtree = indtree_fit
  ),
  metadata = list(
    dmm_mixture_weight = dmm_mixture_weight,
    tree_depth = tree_depth
  )
)

saveRDS(comparison_bundle, out_bundle)
saveRDS(dmm_fit, out_dmm_fit)
saveRDS(indtree_fit, out_indtree_fit)
saveRDS(dmm_posterior, out_dmm_posterior)
saveRDS(indtree_pi_trace, out_indtree_pi)
utils::write.csv(cluster_assignments, out_cluster_csv, row.names = FALSE)
utils::write.csv(timing_df, out_timing_csv, row.names = FALSE)
utils::write.csv(method_summary, out_method_summary_csv, row.names = FALSE)
utils::write.csv(crosstab_vs_cortree, out_crosstab_csv, row.names = FALSE)
utils::write.csv(covariate_assoc_df, out_cov_assoc_csv, row.names = FALSE)
saveRDS(covariate_assoc_list, out_cov_assoc_rds)

grDevices::pdf(out_pdf, width = 14, height = 5)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(
  umap_embedding[, 1],
  umap_embedding[, 2],
  col = cluster_cortree + 1L,
  pch = 19,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = "CorTree"
)
plot(
  umap_embedding[, 1],
  umap_embedding[, 2],
  col = cluster_dmm + 1L,
  pch = 19,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = sprintf("DMM (ARI %.3f)", method_summary$ari_vs_cortree[method_summary$method == "DMM"])
)
plot(
  umap_embedding[, 1],
  umap_embedding[, 2],
  col = cluster_indtree + 1L,
  pch = 19,
  xlab = "UMAP1",
  ylab = "UMAP2",
  main = sprintf("IndTree (ARI %.3f)", method_summary$ari_vs_cortree[method_summary$method == "IndTree"])
)
grDevices::dev.off()

cat("AGP comparison complete\n")
cat("Loaded CorTree bundle:", cortree_bundle_path, "\n")
cat("Saved comparison bundle:", out_bundle, "\n")
cat("Saved cluster assignments:", out_cluster_csv, "\n")
cat("Saved timing summary:", out_timing_csv, "\n")
cat("Saved method summary:", out_method_summary_csv, "\n")
cat("Saved cross-tab summary:", out_crosstab_csv, "\n")
cat("Saved covariate association CSV:", out_cov_assoc_csv, "\n")
cat("Saved covariate association RDS:", out_cov_assoc_rds, "\n")
cat("Saved DMM fit:", out_dmm_fit, "\n")
cat("Saved IndTree fit:", out_indtree_fit, "\n")
cat("Saved DMM posterior:", out_dmm_posterior, "\n")
cat("Saved IndTree pi trace:", out_indtree_pi, "\n")
cat("Saved UMAP comparison PDF:", out_pdf, "\n")
print(timing_df)
print(method_summary)
