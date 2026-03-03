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
  cov_interval = cov_interval
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


# run ind-tree ----------------------------------------------------------

# removed err_precision to prevent singularity
tree_depth <- 6
cutoff_layer <- 3    # Layer for correlated nodes cutoff

burnin <- 100        # Burn-in period
total_iter <- burnin + 50   # Total iterations
c_sigma2_vec <- 10  # Hyperparameter for sigma^2 vector. 1/c is the mean
sigma_mu2 <- 0.1     # Hyperparameter for sigma_mu^2
# Call the CorTree_sampler function

indtree <- CorTree_sampler(X, 
                            # init_Z = Z_pam-1,
                            init_Z = X_rowsum-1,
                            n_clus = n_clus+1, 
                            tree_depth, 
                            cutoff_layer, 
                            total_iter, 
                            burnin, 
                            c_sigma2_vec, 
                            sigma_mu2,
                            warm_start = 0,
                            all_ind = T)
indtree$elapsed
par(mfrow = c(1,1))
plot(indtree$mcmc$loglik,type="l")

# clustering pattern
Z_indtree = apply(indtree$mcmc$Z[,-c(total_iter-burnin)],1,mean);table(Z_indtree)
Z_indtree = round(Z_indtree); table(Z_indtree)

mclust::adjustedRandIndex(Z_true, Z_indtree)

