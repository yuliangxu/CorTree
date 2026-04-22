rm(list = ls())

data_path = "/hpc/group/mastatlab/yx306/DNA_data"
out_dir = "/cwork/yx306/CorTree"
project_dir <- normalizePath(".")
source(file.path(project_dir, "R", "microbiome_help.R"))
source(file.path(project_dir, "R", "DCC_help.R"))
devtools::load_all()

args <- commandArgs(trailingOnly = TRUE)
valid_init_modes <- c("rowsum_quantile", "pwm_quantile", "pam", "tss_dist_quantile")
task_id <- if (length(args) >= 1L) {
  suppressWarnings(as.integer(args[[1L]]))
} else {
  suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1")))
}
if (is.na(task_id) || task_id < 1L || task_id > length(valid_init_modes)) {
  stop("Task id must be in 1:", length(valid_init_modes))
}
init_mode <- valid_init_modes[[task_id]]

file_name <- "NRF1.K562.DNase.counts.mat.rds"
label_name <- "NRF1.K562.sites.chip.labels.tss.dist.rds"
data_dir <- resolve_data_dir(data_path, c(file_name, label_name))

dnase_count_matrix <- readRDS(file.path(data_dir, file_name))
sites_chip_labels <- readRDS(file.path(data_dir, label_name))
dnase_count_matrix <- dnase_count_matrix[, 1:(ncol(dnase_count_matrix) / 2)] +
  dnase_count_matrix[, (ncol(dnase_count_matrix) / 2 + 1):ncol(dnase_count_matrix)]

threshold <- 50
pwm_score_cutoff <- 13
idx_include <- which(
  rowSums(dnase_count_matrix) >= threshold &
    sites_chip_labels$pwm.score >= pwm_score_cutoff &
    sites_chip_labels$strand == "+"
)

X <- as.matrix(dnase_count_matrix[idx_include, , drop = FALSE])
chip_labels <- as.numeric(sites_chip_labels$chip_label[idx_include])
pwm_score <- sites_chip_labels$pwm.score[idx_include]
tss_dist <- sites_chip_labels$tss.dist[idx_include]
row_sum_count <- rowSums(X)
chip_data <- sites_chip_labels$chip[idx_include]

n_clus_bin <- 2L
data_matrix <- scale(X)
t_kmeans <- system.time({
  kmeans_result <- kmeans(data_matrix, centers = n_clus_bin, nstart = 25)
})
Z_kmeans <- as.integer(kmeans_result$cluster)

# DMM ------------------------------------------------
dmm_time <- system.time({
  dmm_fit <- DirichletMultinomial::dmn(X, k = n_clus_bin)
})
dmm_posterior <- extract_dmm_posterior(dmm_fit)
if (nrow(dmm_posterior) != nrow(X)) {
  stop("DMM posterior row count does not match the number of samples.")
}
colnames(dmm_posterior) <- paste0("cluster_", seq_len(ncol(dmm_posterior)) - 1L)
rownames(dmm_posterior) <- rownames(X)
cluster_dmm <- max.col(dmm_posterior, ties.method = "first") - 1L
dmm_mixture_weight <- tryCatch(
  as.numeric(DirichletMultinomial::mixturewt(dmm_fit)),
  error = function(e) rep(NA_real_, n_clus_bin)
)

saveRDS(
  list(
    fit = dmm_fit,
    posterior = dmm_posterior,
    cluster_assignments = cluster_dmm,
    mixture_weight = dmm_mixture_weight,
    elapsed_time = unname(dmm_time[["elapsed"]])
  ),
  file.path(out_dir, paste0("NRF1", "_dmm_results.rds"))
)
# end of DMM ------------------------------------------------

n_clus <- 5L
t_rowsum_init <- system.time({
  init_rowsum <- as.integer(quantile_groups(row_sum_count, k = n_clus) - 1L)
})
t_pwm_init <- system.time({
  init_pwm <- as.integer(quantile_groups(pwm_score, k = n_clus) - 1L)
})
t_tss_dist_init <- system.time({
  init_tss_dist <- as.integer(quantile_groups(tss_dist, k = n_clus) - 1L)
})

library(cluster)
t_pam <- system.time({
  pam_fit <- pam(X, n_clus_bin)
})
Z_pam <- as.integer(pam_fit$clustering)
t_pam_init <- system.time({
  pam_init_fit <- pam(X, n_clus)
})
Z_pam_init <- as.integer(pam_init_fit$clustering - 1L)

library(CENTIPEDE)
t_centipede <- system.time({
  centFit <- fitCentipede(
    Xlist = list(DNase = X),
    Y = cbind(1, pwm_score)
  )
})
Z_centipede <- ifelse(centFit$PostPr > 0.5, 1L, 0L)

init_lookup <- list(
  rowsum_quantile = init_rowsum,
  pwm_quantile = init_pwm,
  pam = Z_pam_init,
  tss_dist_quantile = init_tss_dist
)
init_Z <- init_lookup[[init_mode]]

tree_depth <- 9L
cutoff_layer <- 4L
warm_start <- 0L
burnin <- warm_start + 100L
total_iter <- burnin + 50L
c_sigma2_vec <- 10
sigma_mu2 <- 0.1
cov_interval <- 3L

set.seed(2025)
fit <- CorTree_sampler(
  X = X,
  init_Z = init_Z,
  n_clus = n_clus,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  warm_start = warm_start,
  cov_interval = cov_interval,
  all_ind = FALSE
)

set.seed(2025)
fit_indtree <- CorTree_sampler(
  X = X,
  init_Z = init_Z,
  n_clus = n_clus,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  warm_start = warm_start,
  cov_interval = cov_interval,
  all_ind = TRUE
)

Z_cortree <- hard_cluster_from_chain(fit$mcmc$Z, n_clus = n_clus, burnin = burnin)
Z_indtree <- hard_cluster_from_chain(fit_indtree$mcmc$Z, n_clus = n_clus, burnin = burnin)
pi_trace <- extract_pi_trace(fit, burnin = burnin)
pi_trace_indtree <- extract_pi_trace(fit_indtree, burnin = burnin)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file_tag <- sprintf("NRF1_task%02d_%s", task_id, init_mode)

init_label <- paste0("Z_init (", init_mode, ")")
summary_table <- data.frame(
  Method = c("K-means", "PAM", "CENTIPEDE", init_label, "IndTree", "CorTree"),
  ARI = c(
    adjusted_rand_index(chip_labels, Z_kmeans),
    adjusted_rand_index(chip_labels, Z_pam),
    adjusted_rand_index(chip_labels, Z_centipede),
    adjusted_rand_index(chip_labels, init_Z),
    adjusted_rand_index(chip_labels, Z_indtree),
    adjusted_rand_index(chip_labels, Z_cortree)
  )
)
summary_table <- summary_table[order(-summary_table$ARI), , drop = FALSE]
utils::write.csv(summary_table, file.path(out_dir, paste0(file_tag, "_ari_table.csv")), row.names = FALSE)
writeLines(
  capture.output(print(knitr::kable(summary_table, digits = 3))),
  con = file.path(out_dir, paste0(file_tag, "_ari_summary.txt"))
)

library(patchwork)

max_log_count <- 800
p0 <- hist_by_group_gg(chip_data, chip_labels, title = "chip_labels", y_max = max_log_count)
p1 <- hist_by_group_gg(chip_data, Z_kmeans, title = "Kmeans", y_max = max_log_count)
p2 <- hist_by_group_gg(chip_data, Z_pam, title = "PAM", y_max = max_log_count)
p3 <- hist_by_group_gg(chip_data, Z_centipede, title = "CENTIPEDE", y_max = max_log_count)
p4 <- hist_by_group_gg(chip_data, init_Z, title = init_label, y_max = max_log_count)
p5 <- hist_by_group_gg(chip_data, Z_indtree, title = "IndTree", y_max = max_log_count)
p6 <- hist_by_group_gg(chip_data, Z_cortree, title = "CorTree", y_max = max_log_count)

y_upper_lim <- NULL
m0 <- cluster_mean_gg(X, chip_labels, title = "chip_labels", y_max = y_upper_lim)
m1 <- cluster_mean_gg(X, Z_kmeans, title = "Kmeans", y_max = y_upper_lim)
m2 <- cluster_mean_gg(X, Z_pam, title = "PAM", y_max = y_upper_lim)
m3 <- cluster_mean_gg(X, Z_centipede, title = "CENTIPEDE", y_max = y_upper_lim)
m4 <- cluster_mean_gg(X, init_Z, title = init_label, y_max = y_upper_lim)
m5 <- cluster_mean_gg(X, Z_indtree, title = "IndTree", y_max = y_upper_lim)
m6 <- cluster_mean_gg(X, Z_cortree, title = "CorTree", y_max = y_upper_lim)

combined <- wrap_plots(
  m0, m1, m2, m3, m4, m5, m6,
  p0, p1, p2, p3, p4, p5, p6,
  ncol = 7
)
ggplot2::ggsave(
  filename = file.path(out_dir, paste0(file_tag, "_combined.pdf")),
  plot = combined,
  width = 28,
  height = 8
)
saveRDS(
  list(
    combined = combined,
    cluster_mean = list(
      chip_labels = m0,
      kmeans = m1,
      pam = m2,
      centipede = m3,
      init = m4,
      indtree = m5,
      cortree = m6
    ),
    histogram = list(
      chip_labels = p0,
      kmeans = p1,
      pam = p2,
      centipede = p3,
      init = p4,
      indtree = p5,
      cortree = p6
    )
  ),
  file.path(out_dir, paste0(file_tag, "_plot_artifacts.rds"))
)

result_bundle <- list(
  init_mode = init_mode,
  input = list(
    data_dir = data_dir,
    file_name = file_name,
    label_name = label_name,
    threshold = threshold,
    pwm_score_cutoff = pwm_score_cutoff,
    idx_include = idx_include,
    sampler_args = list(
      n_clus = n_clus,
      tree_depth = tree_depth,
      cutoff_layer = cutoff_layer,
      total_iter = total_iter,
      burnin = burnin,
      warm_start = warm_start,
      c_sigma2_vec = c_sigma2_vec,
      sigma_mu2 = sigma_mu2,
      cov_interval = cov_interval,
      all_ind = FALSE
    )
  ),
  data = list(
    X = X,
    chip_labels = chip_labels,
    pwm_score = pwm_score,
    tss_dist = tss_dist,
    row_sum_count = row_sum_count,
    chip_data = chip_data
  ),
  init = list(
    rowsum_quantile = init_rowsum,
    pwm_quantile = init_pwm,
    pam = Z_pam_init,
    tss_dist_quantile = init_tss_dist,
    selected = init_Z
  ),
  cluster_assignments = list(
    kmeans = Z_kmeans,
    pam = Z_pam,
    centipede = Z_centipede,
    init = init_Z,
    indtree = Z_indtree,
    cortree = Z_cortree
  ),
  fit = list(
    cortree = fit,
    indtree = fit_indtree
  ),
  timing = list(
    kmeans = unname(t_kmeans[["elapsed"]]),
    pam = unname(t_pam[["elapsed"]]),
    centipede = unname(t_centipede[["elapsed"]]),
    init = list(
      rowsum_quantile = unname(t_rowsum_init[["elapsed"]]),
      pwm_quantile = unname(t_pwm_init[["elapsed"]]),
      pam = unname(t_pam_init[["elapsed"]]),
      tss_dist_quantile = unname(t_tss_dist_init[["elapsed"]])
    ),
    indtree = fit_indtree$elapsed,
    cortree = fit$elapsed
  ),
  trace = list(
    cortree_pi_postburnin = pi_trace,
    indtree_pi_postburnin = pi_trace_indtree
  ),
  ari_table = summary_table
)
saveRDS(result_bundle, file.path(out_dir, paste0(file_tag, "_result_bundle.rds")))
saveRDS(pi_trace, file.path(out_dir, paste0(file_tag, "_cortree_pi_trace_postburnin.rds")))
saveRDS(pi_trace_indtree, file.path(out_dir, paste0(file_tag, "_indtree_pi_trace_postburnin.rds")))

cat("NRF1 init comparison complete\n")
cat("Task id:", task_id, "\n")
cat("Init mode:", init_mode, "\n")
cat("Saved task-tagged ARI table:", file.path(out_dir, paste0(file_tag, "_ari_table.csv")), "\n")
cat("Saved task-tagged ARI summary:", file.path(out_dir, paste0(file_tag, "_ari_summary.txt")), "\n")
cat("Saved task-tagged combined PDF:", file.path(out_dir, paste0(file_tag, "_combined.pdf")), "\n")
cat("Saved task-tagged plot artifacts:", file.path(out_dir, paste0(file_tag, "_plot_artifacts.rds")), "\n")
cat("Saved task-tagged result bundle:", file.path(out_dir, paste0(file_tag, "_result_bundle.rds")), "\n")
cat("Saved task-tagged post-burnin pi trace:", file.path(out_dir, paste0(file_tag, "_cortree_pi_trace_postburnin.rds")), "\n")
cat("Saved task-tagged post-burnin pi trace:", file.path(out_dir, paste0(file_tag, "_indtree_pi_trace_postburnin.rds")), "\n")
print(summary_table)
