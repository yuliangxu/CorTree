rm(list = ls())

data_path = "/hpc/group/mastatlab/yx306/DNA_data"
out_dir = "/cwork/yx306/CorTree"
project_dir <- normalizePath(".")
source(file.path(project_dir, "R", "microbiome_help.R"))
source(file.path(project_dir, "R", "DCC_help.R"))
devtools::load_all()

args <- commandArgs(trailingOnly = TRUE)
valid_sensitivity_modes <- c("cutoff_layer_3", "c_sigma2_vec_1", "sigma_mu2_1")
n_tasks <- length(valid_sensitivity_modes)

task_id <- if (length(args) >= 1L) {
  suppressWarnings(as.integer(args[[1L]]))
} else {
  suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1")))
}
if (is.na(task_id) || task_id < 1L || task_id > n_tasks) {
  stop("Task id must be in 1:", n_tasks)
}

init_mode <- "pwm_quantile"
sensitivity_mode <- valid_sensitivity_modes[[task_id]]

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
t_pam_init <- system.time({
  pam_init_fit <- pam(X, n_clus)
})
Z_pam_init <- as.integer(pam_init_fit$clustering - 1L)

init_lookup <- list(
  rowsum_quantile = init_rowsum,
  pwm_quantile = init_pwm,
  pam = Z_pam_init,
  tss_dist_quantile = init_tss_dist
)
init_timing_lookup <- list(
  rowsum_quantile = unname(t_rowsum_init[["elapsed"]]),
  pwm_quantile = unname(t_pwm_init[["elapsed"]]),
  pam = unname(t_pam_init[["elapsed"]]),
  tss_dist_quantile = unname(t_tss_dist_init[["elapsed"]])
)
init_Z <- init_lookup[[init_mode]]

baseline_sampler_args <- list(
  n_clus = n_clus,
  tree_depth = 9L,
  cutoff_layer = 4L,
  warm_start = 0L,
  burnin = 100L,
  total_iter = 150L,
  c_sigma2_vec = 10,
  sigma_mu2 = 0.1,
  cov_interval = 3L,
  all_ind = FALSE
)

sampler_args <- baseline_sampler_args
if (identical(sensitivity_mode, "cutoff_layer_3")) {
  sampler_args$cutoff_layer <- 3L
} else if (identical(sensitivity_mode, "c_sigma2_vec_1")) {
  sampler_args$c_sigma2_vec <- 1
} else if (identical(sensitivity_mode, "sigma_mu2_1")) {
  sampler_args$sigma_mu2 <- 1
} else {
  stop("Unsupported sensitivity mode: ", sensitivity_mode)
}

sensitivity_label <- switch(
  sensitivity_mode,
  cutoff_layer_3 = "cutoff_layer=4",
  c_sigma2_vec_1 = "c_sigma2_vec=1",
  sigma_mu2_1 = "sigma_mu2=1"
)

set.seed(2025)
fit <- CorTree_sampler(
  X = X,
  init_Z = init_Z,
  n_clus = sampler_args$n_clus,
  tree_depth = sampler_args$tree_depth,
  cutoff_layer = sampler_args$cutoff_layer,
  total_iter = sampler_args$total_iter,
  burnin = sampler_args$burnin,
  c_sigma2_vec = sampler_args$c_sigma2_vec,
  sigma_mu2 = sampler_args$sigma_mu2,
  warm_start = sampler_args$warm_start,
  cov_interval = sampler_args$cov_interval,
  all_ind = sampler_args$all_ind
)

Z_cortree <- hard_cluster_from_chain(
  fit$mcmc$Z,
  n_clus = sampler_args$n_clus,
  burnin = sampler_args$burnin
)
pi_trace <- extract_pi_trace(fit, burnin = sampler_args$burnin)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
file_tag <- sprintf(
  "NRF1_sensi_task%02d_%s",
  task_id,
  sensitivity_mode
)

init_label <- paste0("Z_init (", init_mode, ")")
summary_table <- data.frame(
  Method = c(init_label, "CorTree"),
  ARI = c(
    adjusted_rand_index(chip_labels, init_Z),
    adjusted_rand_index(chip_labels, Z_cortree)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  summary_table,
  file.path(out_dir, paste0(file_tag, "_ari_table.csv")),
  row.names = FALSE
)
writeLines(
  c(
    paste("Dataset:", "NRF1"),
    paste("Init mode:", init_mode),
    paste("Sensitivity:", sensitivity_label),
    "",
    capture.output(print(knitr::kable(summary_table, digits = 3)))
  ),
  con = file.path(out_dir, paste0(file_tag, "_ari_summary.txt"))
)

result_bundle <- list(
  init_mode = init_mode,
  sensitivity_mode = sensitivity_mode,
  sensitivity_label = sensitivity_label,
  input = list(
    data_dir = data_dir,
    file_name = file_name,
    label_name = label_name,
    threshold = threshold,
    pwm_score_cutoff = pwm_score_cutoff,
    idx_include = idx_include,
    baseline_sampler_args = baseline_sampler_args,
    sampler_args = sampler_args
  ),
  data = list(
    chip_labels = chip_labels,
    pwm_score = pwm_score,
    tss_dist = tss_dist,
    row_sum_count = row_sum_count
  ),
  init = list(
    rowsum_quantile = init_rowsum,
    pwm_quantile = init_pwm,
    pam = Z_pam_init,
    tss_dist_quantile = init_tss_dist,
    selected = init_Z
  ),
  cluster_assignments = list(
    init = init_Z,
    cortree = Z_cortree
  ),
  fit = list(
    cortree = fit
  ),
  timing = list(
    init = init_timing_lookup,
    cortree = fit$elapsed
  ),
  trace = list(
    cortree_pi_postburnin = pi_trace
  ),
  ari_table = summary_table
)

saveRDS(result_bundle, file.path(out_dir, paste0(file_tag, "_result_bundle.rds")))
saveRDS(
  pi_trace,
  file.path(out_dir, paste0(file_tag, "_cortree_pi_trace_postburnin.rds"))
)

cat("NRF1 sensitivity analysis complete\n")
cat("Task id:", task_id, "\n")
cat("Init mode:", init_mode, "\n")
cat("Sensitivity:", sensitivity_label, "\n")
cat(
  "Saved task-tagged ARI table:",
  file.path(out_dir, paste0(file_tag, "_ari_table.csv")),
  "\n"
)
cat(
  "Saved task-tagged ARI summary:",
  file.path(out_dir, paste0(file_tag, "_ari_summary.txt")),
  "\n"
)
cat(
  "Saved task-tagged result bundle:",
  file.path(out_dir, paste0(file_tag, "_result_bundle.rds")),
  "\n"
)
cat(
  "Saved task-tagged post-burnin pi trace:",
  file.path(out_dir, paste0(file_tag, "_cortree_pi_trace_postburnin.rds")),
  "\n"
)
print(summary_table)
