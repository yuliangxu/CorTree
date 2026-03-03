# CorTree test using the same simulation setup as Sim1_clus.R

devtools::load_all()
library(dplyr)

set.seed(2025+2)

# generate samples from mixtures -----------------------------------------------

n_sample = 400
n_leaf = 1000
n_clus = 2
true_params = list()
true_params$pi = c(0.6,0.4)
Z = sample(1:n_clus, n_sample, prob = true_params$pi, replace = TRUE)
Z_true = Z
X = matrix(0, nrow = n_sample, ncol = n_leaf)
Z_tab = table(Z)

gen_X = function(n, Z){
  weight = rbeta(1, 10, 10)
  n1 = floor(n * weight)
  n2 = n - n1
  switch(Z,
         c(rbeta(n1, 2, 6), rbeta(n2, 6, 2)),
         c(rbeta(n1, 1, 1), rbeta(n2, 3, 3))
  )
}

for(i in 1:n_clus){
  X_list <- replicate(
    n = as.numeric(Z_tab[i]),
    gen_X(sample(4e3:1e4, 1), i)
  )
  hist_breaks = seq(0, 1, length.out = n_leaf + 1)
  counts_list = lapply(X_list, function(x) hist(x, breaks = hist_breaks, plot = FALSE)$counts)
  counts_matrix <- do.call(cbind, counts_list)
  X[which(Z == i), ] = t(counts_matrix)
}

# optional initial clustering --------------------------------------------------

library(cluster)
pamx <- pam(X, 3)
Z_pam = c(pamx$clustering)

# run CorTree -----------------------------------------------------------------

tree_depth <- 6
cutoff_layer <- 4
burnin <- 200
total_iter <- burnin + 50
c_sigma2_vec <- 10
sigma_mu2 <- 0.1
cov_interval <- 5
X_rowsum_raw = rowSums(X)
X_rowsum = dplyr::ntile(X_rowsum_raw, n_clus + 1)
devtools::load_all()
cortree <- CorTree_sampler(X = X,
  n_clus = n_clus + 1L,
  init_Z = X_rowsum - 1,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  warm_start = 0,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  all_ind = FALSE,
  cov_interval = cov_interval,
  save_phi_trace = FALSE,
  save_cluster_cor_trace = FALSE
)

cortree$elapsed

Z_cortree = apply(cortree$mcmc$Z[, -c(total_iter - burnin)], 1, mean)
Z_cortree = round(Z_cortree)
table(Z_cortree)
mclust::adjustedRandIndex(Z_true, Z_cortree)

plot(cortree$mcmc$loglik, type = "l")

# trace plot of pi
pi_trace = t(cortree$mcmc$pi)
matplot(pi_trace, type = "l", lty = 1, col = 1:ncol(pi_trace))
legend("topright", legend = paste0("Cluster ", 1:ncol(pi_trace)), col = 1:ncol(pi_trace), lty = 1)


# cortree with random init ----------------------------------------------------------
devtools::load_all()
burnin <- 30
total_iter <- burnin + 30
user_init_Z_list <- list(X_rowsum - 1L)
fit_multi <- CorTree_sampler_randominit(
  X = X,
  n_clus = n_clus + 1L,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  n_start = 5L,
  init_Z_list = user_init_Z_list,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  all_ind = FALSE,
  cov_interval = cov_interval,
  save_phi_trace = FALSE,
  save_cluster_cor_trace = FALSE
)
fit_multi$elapsed_all
Z_cortree_multi = apply(fit_multi$best_fit$mcmc$Z[, -c(total_iter - burnin)], 1, mean)
Z_cortree_multi = round(Z_cortree_multi)
table(Z_cortree_multi)
mclust::adjustedRandIndex(Z_true, Z_cortree_multi)

# pi trace subplot for all starts
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), oma = c(0, 0, 2, 0))
n_user_init <- length(user_init_Z_list)
best_start_id <- fit_multi$best_start_id
ari_by_start <- rep(NA_real_, length(fit_multi$all_fit))
for (s in seq_along(fit_multi$all_fit)) {
  Z_chain_s <- fit_multi$all_fit[[s]]$mcmc$Z
  Z_hat_s <- round(rowMeans(Z_chain_s))
  ari_by_start[s] <- mclust::adjustedRandIndex(Z_true, Z_hat_s)
}
for (s in seq_along(fit_multi$pi_mcmc_all)) {
  pi_trace_s <- t(fit_multi$pi_mcmc_all[[s]])
  matplot(pi_trace_s, type = "l", lty = 1, col = seq_len(ncol(pi_trace_s)),
          xlab = "Iter", ylab = "pi", main = paste("Start", s))
  subtitle_parts <- character(0)
  subtitle_parts <- c(subtitle_parts, if (s <= n_user_init) "User-given init_Z" else "Random init_Z")
  if (s == best_start_id) subtitle_parts <- c(subtitle_parts, "Best fit chain")
  subtitle_parts <- c(subtitle_parts, sprintf("ARI=%.3f", ari_by_start[s]))
  subtitle_s <- paste(subtitle_parts, collapse = " | ")
  mtext(subtitle_s, side = 3, line = 0.1, cex = 0.8)
}
mtext("CorTree pi Trace by Initialization", outer = TRUE, cex = 1.1)

# loglik trace subplot for all starts
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), oma = c(0, 0, 2, 0))
for (s in seq_along(fit_multi$pi_mcmc_all)) {
  loglik_s <- fit_multi$all_fit[[s]]$mcmc$loglik
  plot(loglik_s, type = "l", col = "black",
       xlab = "Iter", ylab = "loglik", main = paste("Start", s))
  abline(v = burnin, col = "red", lty = 2)
  subtitle_parts <- character(0)
  subtitle_parts <- c(subtitle_parts, if (s <= n_user_init) "User-given init_Z" else "Random init_Z")
  if (s == best_start_id) subtitle_parts <- c(subtitle_parts, "Best fit chain")
  subtitle_s <- paste(subtitle_parts, collapse = " | ")
  mtext(subtitle_s, side = 3, line = 0.1, cex = 0.8)
}
mtext("CorTree loglik Trace by Initialization", outer = TRUE, cex = 1.1)

# held-out log predictive density model selection ------------------------------------
devtools::load_all()
set.seed(2025 + 200)
train_idx <- sample(seq_len(nrow(X)), size = floor(0.8 * nrow(X)), replace = FALSE)
test_idx <- setdiff(seq_len(nrow(X)), train_idx)

burnin <- 50
total_iter <- burnin + 50
X_rowsum_train <- dplyr::ntile(rowSums(X[train_idx, , drop = FALSE]), n_clus + 1L)
heldout_fit <- CorTree_sampler_randominit_heldout(
  X = X,
  train_idx = train_idx,
  test_idx = test_idx,
  n_clus = n_clus + 1L,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  n_start = 5L,
  init_Z_list = list(X_rowsum_train - 1L),
  warm_start = 0L,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  all_ind = FALSE,
  cov_interval = cov_interval,
  save_phi_trace = FALSE,
  save_cluster_cor_trace = FALSE,
  n_phi_mc = 30L
)

heldout_fit$best_start_id_heldout
heldout_fit$heldout_hlpd_by_start

Z_train_hat <- round(rowMeans(heldout_fit$best_fit_heldout$mcmc$Z))
mclust::adjustedRandIndex(Z_true[train_idx], Z_train_hat)

# diagnostics for held-out model selection
old_par_heldout <- par(no.readonly = TRUE)
n_user_init_heldout <- 1L
best_start_id_heldout <- heldout_fit$best_start_id_heldout
hlpd_by_start <- heldout_fit$heldout_hlpd_by_start
ari_by_start_heldout <- rep(NA_real_, length(heldout_fit$fit_multi$all_fit))
for (s in seq_along(heldout_fit$fit_multi$all_fit)) {
  Z_chain_s <- heldout_fit$fit_multi$all_fit[[s]]$mcmc$Z
  Z_hat_s <- round(rowMeans(Z_chain_s))
  ari_by_start_heldout[s] <- mclust::adjustedRandIndex(Z_true[train_idx], Z_hat_s)
}

# (1) pi trace for each random init with held-out HLPD in subtitle
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), oma = c(0, 0, 2, 0))
for (s in seq_along(heldout_fit$fit_multi$pi_mcmc_all)) {
  pi_trace_s <- t(heldout_fit$fit_multi$pi_mcmc_all[[s]])
  matplot(pi_trace_s, type = "l", lty = 1, col = seq_len(ncol(pi_trace_s)),
          xlab = "Iter", ylab = "pi", main = paste("Held-out Start", s))
  subtitle_parts <- character(0)
  subtitle_parts <- c(subtitle_parts, if (s <= n_user_init_heldout) "User-given init_Z" else "Random init_Z")
  if (s == best_start_id_heldout) subtitle_parts <- c(subtitle_parts, "Best held-out fit")
  subtitle_parts <- c(subtitle_parts, sprintf("HLPD=%.2f", hlpd_by_start[s]))
  subtitle_parts <- c(subtitle_parts, sprintf("ARI=%.3f", ari_by_start_heldout[s]))
  mtext(paste(subtitle_parts, collapse = " | "), side = 3, line = 0.1, cex = 0.8)
}
mtext("Held-out Selection: pi Trace by Initialization", outer = TRUE, cex = 1.1)

# (2) loglik trace for each random init with held-out HLPD in subtitle
par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), oma = c(0, 0, 2, 0))
for (s in seq_along(heldout_fit$fit_multi$all_fit)) {
  loglik_s <- heldout_fit$fit_multi$all_fit[[s]]$mcmc$loglik
  plot(loglik_s, type = "l", col = "black",
       xlab = "Iter", ylab = "loglik", main = paste("Held-out Start", s))
  abline(v = burnin, col = "red", lty = 2)
  subtitle_parts <- character(0)
  subtitle_parts <- c(subtitle_parts, if (s <= n_user_init_heldout) "User-given init_Z" else "Random init_Z")
  if (s == best_start_id_heldout) subtitle_parts <- c(subtitle_parts, "Best held-out fit")
  subtitle_parts <- c(subtitle_parts, sprintf("HLPD=%.2f", hlpd_by_start[s]))
  subtitle_parts <- c(subtitle_parts, sprintf("ARI=%.3f", ari_by_start_heldout[s]))
  mtext(paste(subtitle_parts, collapse = " | "), side = 3, line = 0.1, cex = 0.8)
}
mtext("Held-out Selection: loglik Trace by Initialization", outer = TRUE, cex = 1.1)
par(old_par_heldout)
