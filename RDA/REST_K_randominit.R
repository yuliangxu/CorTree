
# data_path = "/hpc/group/mastatlab/yx306/DNA_data/"
data_path = '/Users/yuliangxu/Library/Mobile Documents/com~apple~CloudDocs/MyDrive/Research/Tree/code/AutoTree/data/'
list.files(path = data_path)

source("./R/microbiome_help.R")
file_name_K_new = "REST.K562.DNase.counts.mat.rds"
label_K_new = "REST.K562.sites.chip.labels.rds"



dnase_count_matrix <- readRDS(paste0(data_path,file_name_K_new))
sites_chip_labels <- readRDS(paste0(data_path,label_K_new))
dnase_count_matrix <- dnase_count_matrix[,1:(ncol(dnase_count_matrix)/2)] + dnase_count_matrix[,(ncol(dnase_count_matrix)/2+1):ncol(dnase_count_matrix)]
pos_idx <- which(sites_chip_labels$chip_label == 1)
neg_idx <- which(sites_chip_labels$chip_label == 0)




threshold = 100
pwm_score_cutoff = 13


idx_include = which(rowSums(dnase_count_matrix) >= threshold & sites_chip_labels$pwm.score >= pwm_score_cutoff )
# idx_include = which(rowSums(dnase_count_matrix) >= threshold & sites_chip_labels$pwm.score >= pwm_score_cutoff & 
#                       sites_chip_labels$chr %in% c("chr1","chr2","chr3","chr4","chr5","chr6"))
length(idx_include)

set.seed(2025)
all_include = idx_include

X = as.matrix(dnase_count_matrix[all_include,]); dim(X)

chip_labels = as.numeric(sites_chip_labels$chip_label[all_include])
pos_idx = which(chip_labels ==1)
neg_idx = which(chip_labels ==0)

plot_pos_neg_profiles(X, chip_labels,
                      title = "Aggregate DNase-seq profiles around REST motifs in K562")
summary(rowSums(X))

n_clus = 5
y_upper_lim = 15
in_tree_depth = 9

pwm_score = sites_chip_labels$pwm.score[all_include]

quantile_groups <- function(x, k) {
  qs <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    labels = seq_len(k)
  )
}
# initial_value = 1*(pwm_score >= median(pwm_score) )
# initial_value = as.numeric(quantile_groups(pwm_score,k=n_clus)) # PWM-based initial value
total_count = rowSums(X)
initial_value = as.numeric(quantile_groups(total_count,k=n_clus)) # TotalCount-based initial value

table(initial_value,chip_labels)
mclust::adjustedRandIndex(initial_value,chip_labels)

n_clus_bin = 2

# Perform k-means clustering ------------------------------------------------




data_matrix <- scale(X)
kmeans_result <- kmeans(data_matrix, centers = n_clus_bin, nstart = 25)
Z_kmeans = kmeans_result$cluster;
plot_pos_neg_profiles(X,Z_kmeans,title="kmeans")
mclust::adjustedRandIndex(chip_labels,Z_kmeans)

# PAM -----------------------------------------------------------------------

library(cluster)

pamx <- pam(X, n_clus_bin)
Z_pam = c(pamx$clustering)


plot_pos_neg_profiles(X,Z_pam,title="PAM")

mclust::adjustedRandIndex(chip_labels,Z_pam)



# # run ind-tree --------------------------------------------------------
# 
# # # removed err_precision to prevent singularity
# # 
# tree_depth <- in_tree_depth
# cutoff_layer <- 3; print(paste("num of cor nodes",2^(cutoff_layer+1)-1))    # Layer for correlated nodes cutoff
# warm_start <- 0
# burnin <- warm_start + 50  # Burn-in period
# total_iter <- burnin + 50;print(total_iter)  # Total iterations
# c_sigma2_vec <- 10  # Hyperparameter for sigma^2 vector. 1/c is the mean
# sigma_mu2 <- 0.1     # Hyperparameter for sigma_mu^2
# cov_interval <- 5
# # Call the AutoTree_sampler function
# devtools::load_all()
# set.seed(2025)
# ind_tree <- CorTree_sampler(X,
#                              # init_Z = Z_pam-1,
#                              init_Z = initial_value_pwm-1,
#                              n_clus = n_clus,
#                              tree_depth,
#                              cutoff_layer,
#                              total_iter,
#                              burnin,
#                              c_sigma2_vec,
#                              sigma_mu2,
#                              warm_start = warm_start,
#                              cov_interval = cov_interval,
#                              all_ind = T) # set sn_cor, smallest number in each cluster
# 
# ind_tree$elapsed
# par(mfrow = c(1,1))
# plot(ind_tree$mcmc$loglik[10:total_iter],type="l")
# 
# Z_indtree = apply(ind_tree$mcmc$Z[,-c(total_iter-burnin)],1,mean);table(Z_indtree)
# mclust::adjustedRandIndex(chip_labels,Z_indtree)
# # 
# p1 <- cluster_mean_gg(X,Z_indtree, title = "Z_indtree")
# p2_pwm <- cluster_mean_gg(X,initial_value_pwm-1, title = "PWM (init)")
# p3 <- cluster_mean_gg(X,chip_labels, title = "ChIP")
# 
# combined <- patchwork::wrap_plots(
#   p1, p2_pwm, p3,
#   ncol = 3
# )
# 
# print(combined)


# try CENTIPEDE -----------------------------------------------------------

# install.packages("CENTIPEDE",
#                  repos="http://R-Forge.R-project.org")
library(CENTIPEDE)



centFit <- fitCentipede(
  Xlist = list(DNase = as.matrix(X)),
  Y     = cbind(1, pwm_score)
  # Y     = cbind(1, total_count)
  
  
)

# posterior binding probabilities (i.e. cluster assignments)
post_probs <- centFit$PostPr
cluster_assignments <- ifelse(post_probs > 0.5, 1, 0)
Z_centipede = cluster_assignments

mclust::adjustedRandIndex(chip_labels, Z_centipede) # 0.1577784



# using count as the predictor

centFit_count <- fitCentipede(
  Xlist = list(DNase = as.matrix(X)),
#   Y     = cbind(1, pwm_score)
  Y     = cbind(1, total_count)
  
  
)
post_probs_count <- centFit_count$PostPr
cluster_assignments_count <- ifelse(post_probs_count > 0.5, 1, 0)          
Z_centipede_count = cluster_assignments_count
mclust::adjustedRandIndex(chip_labels, Z_centipede_count) #0.1462061


# run cor-tree ------------------------------------------------------------

n_clus = 5
tree_depth <- in_tree_depth # final result
print(paste("num of leaf nodes",2^(tree_depth-1))) 
cutoff_layer <- 6; print(paste("num of cor nodes",2^(cutoff_layer+1)-1))    # Layer for correlated nodes cutoff
warm_start <- 0
burnin <- warm_start + 20  # Burn-in period
total_iter <- burnin + 20; print(total_iter)  # Total iterations
c_sigma2_vec <- 10  # Hyperparameter for sigma^2 vector. 1/c is the mean, cannot be too large
sigma_mu2 <- 0.1     # Hyperparameter for sigma_mu^2
# cov_interval <- 5 # for data_G only
cov_interval <- 5 # for data_G only
# Call the AutoTree_sampler function
devtools::load_all()

set.seed(2026)
initial_value_count = as.numeric(quantile_groups(total_count,k=n_clus)) 
initial_value_pwm = as.numeric(quantile_groups(pwm_score,k=n_clus)) 

train_idx <- sample(seq_len(nrow(X)), size = floor(0.8 * nrow(X)), replace = FALSE)
test_idx <- setdiff(seq_len(nrow(X)), train_idx)

heldout_fit <- CorTree_sampler_randominit_heldout(
  X = X,
  train_idx = train_idx,
  test_idx = test_idx,
  n_clus = n_clus,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  n_start = 9L,
  init_Z_list = list(count_init = as.integer(initial_value_count - 1L), 
                     pwm_init = as.integer(initial_value_pwm - 1L)),
  warm_start = warm_start,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  all_ind = FALSE,
  cov_interval = cov_interval,
  save_phi_trace = FALSE,
  save_cluster_cor_trace = FALSE,
  n_phi_mc = 20L,
  refit_full_data = TRUE
)

print(paste("best start by held-out HLPD:", heldout_fit$best_start_id_heldout))
print(heldout_fit$heldout_hlpd_by_start)

heldout_fit$train_elapsed_by_start
heldout_fit$heldout_hlpd_elapsed_by_start
heldout_fit$full_fit_elapsed


result <- heldout_fit$full_fit_from_combined_init
if (is.null(result)) {
  stop("`full_fit_from_combined_init` is NULL. Set `refit_full_data = TRUE`.")
}

# Full-data converged labels from post-burnin samples
Z_chain_full <- result$mcmc$Z
n_keep <- total_iter - burnin
if (ncol(Z_chain_full) > n_keep) {
  Z_chain_full <- Z_chain_full[, (burnin + 1L):ncol(Z_chain_full), drop = FALSE]
}
Z_cortree <- as.integer(round(rowMeans(Z_chain_full)))
Z_cortree <- pmax(0L, pmin(as.integer(n_clus) - 1L, Z_cortree))

table(Z_cortree)
mclust::adjustedRandIndex(chip_labels, Z_cortree)

heldout_fit$best_start_id_heldout
heldout_fit$heldout_hlpd_by_start

Z_true = chip_labels
Z_train_hat <- round(rowMeans(heldout_fit$best_fit_heldout$mcmc$Z))
mclust::adjustedRandIndex(Z_true[train_idx], Z_train_hat)

# diagnostics for held-out model selection
old_par_heldout <- par(no.readonly = TRUE)
best_start_id_heldout <- heldout_fit$best_start_id_heldout
hlpd_by_start <- heldout_fit$heldout_hlpd_by_start
n_start_heldout <- length(heldout_fit$fit_multi$all_fit)
init_names_heldout <- names(heldout_fit$fit_multi$init_Z_all)
if (is.null(init_names_heldout)) {
  init_names_heldout <- rep("", n_start_heldout)
}
init_label_by_start <- vapply(
  seq_len(n_start_heldout),
  function(i) {
    nm <- init_names_heldout[i]
    if (!is.na(nm) && nzchar(nm)) {
      return(sprintf("Init: %s", nm))
    }
    "Init: random"
  },
  character(1)
)
ari_by_start_heldout <- rep(NA_real_, length(heldout_fit$fit_multi$all_fit))
for (s in seq_along(heldout_fit$fit_multi$all_fit)) {
  Z_chain_s <- heldout_fit$fit_multi$all_fit[[s]]$mcmc$Z
  Z_hat_s <- round(rowMeans(Z_chain_s))
  ari_by_start_heldout[s] <- mclust::adjustedRandIndex(Z_true[train_idx], Z_hat_s)
}

# (1) pi trace for each random init with held-out HLPD in subtitle
par(mfrow = c(3, 3), mar = c(3, 3, 4.5, 1), oma = c(0, 0, 2, 0))
for (s in seq_along(heldout_fit$fit_multi$pi_mcmc_all)) {
  pi_trace_s <- t(heldout_fit$fit_multi$pi_mcmc_all[[s]])
  matplot(pi_trace_s, type = "l", lty = 1, col = seq_len(ncol(pi_trace_s)),
          xlab = "Iter", ylab = "pi", main = "")
  title(main = paste("Held-out Start", s), line = 2.2)
  subtitle_parts <- character(0)
  subtitle_parts <- c(subtitle_parts, init_label_by_start[s])
  if (s == best_start_id_heldout) subtitle_parts <- c(subtitle_parts, "Best held-out fit")
  subtitle_parts <- c(subtitle_parts, sprintf("HLPD=%.2f", hlpd_by_start[s]))
  mtext(paste(subtitle_parts, collapse = " | "), side = 3, line = 1.0, cex = 0.8)
  mtext(sprintf("ARI=%.3f", ari_by_start_heldout[s]), side = 3, line = 0.2, cex = 0.8)
}
mtext("Held-out Selection: pi Trace by Initialization", outer = TRUE, cex = 1.1)

# (2) loglik trace for each random init with held-out HLPD in subtitle
par(mfrow = c(3, 3), mar = c(3, 3, 4.5, 1), oma = c(0, 0, 2, 0))
for (s in seq_along(heldout_fit$fit_multi$all_fit)) {
  loglik_s <- heldout_fit$fit_multi$all_fit[[s]]$mcmc$loglik
  plot(loglik_s, type = "l", col = "black",
       xlab = "Iter", ylab = "loglik", main = "")
  title(main = paste("Held-out Start", s), line = 2.2)
  abline(v = burnin, col = "red", lty = 2)
  subtitle_parts <- character(0)
  subtitle_parts <- c(subtitle_parts, init_label_by_start[s])
  if (s == best_start_id_heldout) subtitle_parts <- c(subtitle_parts, "Best held-out fit")
  subtitle_parts <- c(subtitle_parts, sprintf("HLPD=%.2f", hlpd_by_start[s]))
  subtitle_parts <- c(subtitle_parts, sprintf("ARI=%.3f", ari_by_start_heldout[s]))
  mtext(paste(subtitle_parts, collapse = " | "), side = 3, line = 0.1, cex = 0.8)
}
mtext("Held-out Selection: loglik Trace by Initialization", outer = TRUE, cex = 1.1)
par(old_par_heldout)


# run cor-tree (PWM init on full data only) ---------------------------------

result_pwm_full <- CorTree_sampler(
  X = X,
  init_Z = as.integer(initial_value_pwm - 1L),
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

result_pwm_full$elapsed
Z_chain_pwm_full <- result_pwm_full$mcmc$Z
if (ncol(Z_chain_pwm_full) > (total_iter - burnin)) {
  Z_chain_pwm_full <- Z_chain_pwm_full[, (burnin + 1L):ncol(Z_chain_pwm_full), drop = FALSE]
}
Z_cortree_pwm_full <- as.integer(round(rowMeans(Z_chain_pwm_full)))
Z_cortree_pwm_full <- pmax(0L, pmin(as.integer(n_clus) - 1L, Z_cortree_pwm_full))
table(Z_cortree_pwm_full)
mclust::adjustedRandIndex(chip_labels, Z_cortree_pwm_full)


# compare with other methods with held-out model selection ------------------------------------------------------------

par(mfrow = c(1,1))
chip_data = sites_chip_labels$chip[all_include]

p1 <- cluster_mean_gg(X,Z_cortree, title = "Z_cortree")
p2_pwm <- cluster_mean_gg(X,initial_value_pwm-1, title = "PWM (init)")
p2_count <- cluster_mean_gg(X,initial_value_count-1, title = "Count (init)")
p3 <- cluster_mean_gg(X,chip_labels, title = "ChIP")
p4 <- cluster_mean_gg(X,Z_centipede, title = "CENTIPEDE")

max_log_count = 300
h1 <- hist_by_group_gg(chip_data,Z_cortree, title = "Cor-tree",y_max = max_log_count)
h2_pwm <- hist_by_group_gg(chip_data,initial_value_pwm-1, title = "PWM (init)",y_max = max_log_count)
h2_count <- hist_by_group_gg(chip_data,initial_value_count-1, title = "Count (init)",y_max = max_log_count)
h3 <- hist_by_group_gg(chip_data,chip_labels, title = "chip_labels",y_max = max_log_count)
h4 <- hist_by_group_gg(chip_data,Z_centipede, title = "CENTIPEDE",y_max = max_log_count)

# combined <- patchwork::wrap_plots(
#   p1, p2_pwm, p3, p4,
#   h1, h2_pwm, h3, h4,
#   ncol = 4
# )
# 
# print(combined)

combined <- patchwork::wrap_plots(
  p1, p2_count, p3, p4,
  h1, h2_count, h3, h4,
  ncol = 4
)

print(combined)



# pi_mcmc = readRDS(file.path(data_path,"cortree_CTCF_K_DNase_sensi_nmcmc30.rds"))
# pi_mcmc = readRDS(file.path(data_path,"cortree_NRF1_K_DNase_sensi_nmcmc30.rds"))

# summarize ---------------------------------------------------------------

Z_cortree3 = readRDS(file.path(data_path,"cortree_REST_K_DNase_sensi_init_kmeans.rds"))
Z_cortree2 = readRDS(file.path(data_path,"cortree_REST_K_DNase_sensi_init_PAM.rds"))
Z_cortree1 = readRDS(file.path(data_path,"cortree_REST_K_DNase_sensi_init_CENTIPEDE.rds"))
Z_cortree = readRDS(file.path(data_path,"cortree_REST_K_DNase.rds"))


# check histogram with Chip data
chip_data = sites_chip_labels$chip[all_include]
source("./R/microbiome_help.R")
library(patchwork)
library(ggplot2)
# every ggplot you build will now use base_size = 14
theme_set(
  theme_minimal(base_size = 14) +
    theme(
      # if you want to tweak axis‐text, legend, etc. independently:
      axis.text    = element_text(size = 12),
      axis.title   = element_text(size = 14),
      plot.title   = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
)

max_log_count = 300
p0 <- hist_by_group_gg(chip_data,chip_labels, title = "chip_labels",y_max = max_log_count)
p1 <- hist_by_group_gg(chip_data,Z_cortree1, title = "Cor-tree (K=6)",y_max = max_log_count)
p2 <- hist_by_group_gg(chip_data,Z_cortree2, title = "Cor-tree (c=1)",y_max = max_log_count)
p3 <- hist_by_group_gg(chip_data,Z_cortree3, title = "Cor-tree (sig2=0.1)",y_max = max_log_count)
p4 <- hist_by_group_gg(chip_data,Z_cortree, title = "Cor-tree",y_max = max_log_count)



# 
# combined <- (p0+p1 )/( p2 + p3)   # or: p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# print(combined)
y_upper_lim = 5
m0 = cluster_mean_gg(X,chip_labels, title = "chip_labels",y_max=y_upper_lim)
m1 <- cluster_mean_gg(X,Z_cortree1, title = "Cor-tree (K=6)",y_max = y_upper_lim)
m2 <- cluster_mean_gg(X,Z_cortree2, title = "Cor-tree (c=1)",y_max = y_upper_lim)
m3 <- cluster_mean_gg(X,Z_cortree3, title = "Cor-tree (sig2=0.1)",y_max = y_upper_lim)
m4 = cluster_mean_gg(X,Z_cortree, title = "Cor-tree",y_max=y_upper_lim)

# combined <- (m0+m1 )/( m2 + m3)   # or: p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# print(combined)

# combined <- (p0 +p1+p2+p3 )/(m0+m1+ m2 + m3)  # or: p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# print(combined)

combined <- wrap_plots(
  m0, m1,
  m2, m3, m4,
  p0, p1,
  p2, p3, p4,
  ncol = 5
)

print(combined)
