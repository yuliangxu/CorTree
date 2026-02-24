# CorTree test using the same simulation setup as Sim1_clus.R

devtools::load_all()
library(dplyr)

set.seed(2025+24)

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
burnin <- 50
total_iter <- burnin + 30
c_sigma2_vec <- 10
sigma_mu2 <- 0.1
cov_interval <- 3
X_rowsum_raw = rowSums(X)
X_rowsum = dplyr::ntile(X_rowsum_raw, n_clus + 1)
devtools::load_all()
cortree <- CorTree_sampler(X = X,
  n_clus = n_clus + 1L,
  init_Z = X_rowsum-1,
  tree_depth = tree_depth,
  cutoff_layer = cutoff_layer,
  total_iter = total_iter,
  burnin = burnin,
  warm_start = 0,
  c_sigma2_vec = c_sigma2_vec,
  sigma_mu2 = sigma_mu2,
  all_ind = FALSE,
  cov_interval = cov_interval
)

cortree$elapsed

# clustering pattern -----------------------------------------------------------

Z_cortree = apply(cortree$mcmc$Z[, -c(total_iter - burnin)], 1, mean)
Z_cortree = round(Z_cortree)
table(Z_cortree)

mclust::adjustedRandIndex(Z_true, Z_cortree)

# check correlated Sigma -------------------------------------------------------

n_cov_mcmc = length(cortree$mcmc$Sigma_inv) - 1
Sigma_inv_mcmc = array(NA, dim = c(n_cov_mcmc, dim(cortree$mcmc$Sigma_inv[[1]])))
for(m in 1:n_cov_mcmc){
  Sigma_inv_mcmc[m, , , ] = cortree$mcmc$Sigma_inv[[m]]
}

Sigma_all_clus = array(NA, dim = c(dim(cortree$mcmc$Sigma_inv[[1]])))
for(k in 1:(n_clus + 1)){
  Sigma_inv_mean = apply(Sigma_inv_mcmc[, , , k], 2:3, mean)
  Sigma_all_clus[, , k] = chol2inv(chol(Sigma_inv_mean))
}

if ((n_clus + 1) >= 3) {
  plot_two_mat(Sigma_all_clus[, , 2], Sigma_all_clus[, , 3])
}

# compare phi mean to empirical phi -------------------------------------------

phi_mean_cor = apply(cortree$mcmc$phi, c(1, 2), mean)
X_tree = construct_tree(X, tree_depth)
empirical_phi = X_tree$empirical_phi

par(mfrow = c(5, 3))
for(j in 1:15){
  plot(phi_mean_cor[, j], empirical_phi[, j], col = Z_true, asp = 1, main = paste0("phi_j", j))
  abline(0, 1)
}
for(j in 16:30){
  plot(phi_mean_cor[, j], empirical_phi[, j], col = Z_true, asp = 1, main = paste0("phi_j", j))
  abline(0, 1)
}
for(j in 31:45){
  plot(phi_mean_cor[, j], empirical_phi[, j], col = Z_true, asp = 1, main = paste0("phi_j", j))
  abline(0, 1)
}
