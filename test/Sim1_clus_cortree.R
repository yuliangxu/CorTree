project_dir <- normalizePath(".")
source(file.path(project_dir, "R", "DCC_help.R"))

out_path = "/cwork/yx306/CorTree"

# Standalone analysis block:
# set run_analysis_only <- TRUE to build the summary table without rerunning simulations.
run_analysis_only <- FALSE

if (run_analysis_only) {
  files <- list.files(
    out_path,
    pattern = "^Sim1_clus_n[0-9]+_rep[0-9]+\\.rds$",
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("No simulation result files found in: ", out_path)
  }

  method_names <- c("K-means", "PAM", "DMM", "Ind-tree", "Cor-tree")
  sample_sizes <- c(200, 400, 600)
  n_method <- length(method_names)

  summary_list <- lapply(sample_sizes, function(n_sample) {
    pattern_n <- paste0("^Sim1_clus_n", n_sample, "_rep[0-9]+\\.rds$")
    files_n <- files[grepl(pattern_n, basename(files))]
    if (length(files_n) == 0L) {
      return(NULL)
    }

    ari_store <- matrix(NA_real_, nrow = length(files_n), ncol = n_method)
    time_store <- matrix(NA_real_, nrow = length(files_n), ncol = n_method)
    keep_idx <- 0L

    for (path in files_n) {
      obj <- readRDS(path)
      ari_mat <- if (!is.null(obj$ari)) as.matrix(obj$ari) else if (!is.null(obj$all_ARI)) as.matrix(obj$all_ARI) else NULL
      time_mat <- if (!is.null(obj$time)) as.matrix(obj$time) else if (!is.null(obj$all_time)) as.matrix(obj$all_time) else NULL
      if (is.null(ari_mat) || is.null(time_mat)) {
        rm(obj, ari_mat, time_mat)
        next
      }
      if (ncol(ari_mat) < n_method || ncol(time_mat) < n_method) {
        rm(obj, ari_mat, time_mat)
        next
      }

      keep_idx <- keep_idx + 1L
      ari_store[keep_idx, ] <- ari_mat[1, seq_len(n_method)]
      time_store[keep_idx, ] <- time_mat[1, seq_len(n_method)]
      rm(obj, ari_mat, time_mat)
    }

    if (keep_idx == 0L) {
      return(NULL)
    }

    make_sim_summary_table(
      ari_mat = ari_store[seq_len(keep_idx), , drop = FALSE],
      time_mat = time_store[seq_len(keep_idx), , drop = FALSE],
      method_names = method_names,
      experiment_label = as.character(n_sample)
    )
  })

  summary_list <- Filter(Negate(is.null), summary_list)
  if (length(summary_list) == 0L) {
    stop("No readable simulation result files were found for summary mode.")
  }

  summary_table <- do.call(rbind, summary_list)
  rownames(summary_table) <- NULL

  summary_csv <- file.path(out_path, "Sim1_clus_summary_table.csv")
  summary_tex <- file.path(out_path, "Sim1_clus_summary_table.tex")
  utils::write.csv(summary_table, summary_csv, row.names = FALSE)
  writeLines(build_sim_summary_latex_table(summary_table), summary_tex)

  print(summary_table)
  cat("Saved summary CSV:", summary_csv, "\n")
  cat("Saved LaTeX table:", summary_tex, "\n")
  quit(save = "no")
}

n_rep = Sys.getenv("SLURM_ARRAY_TASK_ID")
i_rep = as.numeric(n_rep)

devtools::load_all()
# devtools::install()

if (!requireNamespace("DirichletMultinomial", quietly = TRUE)) {
  stop(
    "Sim1_clus_cortree.R needs package 'DirichletMultinomial'. ",
    "Please install it before running the simulation."
  )
}


set_tree_depth = 6

library(dplyr)
set.seed(2025)

# generate samples from mixtures -----------------------------------------------

# n_sample = 600 # 400
if(i_rep <=100){
  n_sample = 200
  
}else if(i_rep <= 200){
  n_sample = 400
  i_rep = i_rep - 100
}else if(i_rep <= 300){
  n_sample = 600
  i_rep = i_rep - 200
}

n_leaf = 1000
n_clus = 2
true_params = list()
true_params$pi = c(0.6,0.4)
Z = sample(1:n_clus, n_sample, prob = true_params$pi, replace=T)
Z_true = Z
# Z = sort(Z)
X = matrix(0, nrow = n_sample, ncol = n_leaf)
Z_tab = table(Z)

gen_X = function(n,Z){
  weight = rbeta(1,10,10)
  n1 = floor(n*weight)
  n2 = n - n1
  switch(Z,
         c(rbeta(n1, 2, 6), rbeta(n2, 6, 2)),
         c(rbeta(n1, 1, 1), rbeta(n2, 3, 3)) 
  )
}

for(i in 1:n_clus){
  X_list <- replicate(
    n = as.numeric(Z_tab[i]),
    gen_X(sample(1e3:5e3, 1),i)
  )
  hist_breaks = seq(0,1,length.out = n_leaf+1)
  counts_list = lapply(X_list, function(x) hist(x, breaks = hist_breaks, plot = FALSE)$counts)
  counts_matrix <- do.call(cbind, counts_list)
  X[which(Z==i),] = t(counts_matrix)
}

X <- matrix(
  as.integer(round(X)),
  nrow = nrow(X),
  ncol = ncol(X),
  dimnames = list(
    paste0("sample_", seq_len(nrow(X))),
    paste0("taxon_", seq_len(ncol(X)))
  )
)

# run knn ---------------------------------------------------------------------

t0 <- Sys.time()
data_matrix <- scale(X)
kmeans_result <- kmeans(data_matrix, centers = 3, nstart = 25)
Z = kmeans_result$cluster; Z_kmeans = Z
table(Z)
clus_label = which(table(Z)>10)
t1 <- Sys.time()
elapsed_knn <- as.numeric(difftime(t1, t0, units = "secs"))



cbind(Z_kmeans, Z_true)[1:10,]
mclust::adjustedRandIndex(Z_true, Z_kmeans)

# run PAM ---------------------------------------------------------------
library(cluster)
t0 <- Sys.time()
pamx <- pam(X, 3)
names(pamx)
length(pamx$clustering)
Z_pam = c(pamx$clustering); table(Z_pam)
t1 <- Sys.time()
elapsed_PAM <- as.numeric(difftime(t1, t0, units = "secs"))

mclust::adjustedRandIndex(Z_true, Z_pam)

# run DMM ---------------------------------------------------------------
t0 <- Sys.time()
dmm_fit <- DirichletMultinomial::dmn(X, k = n_clus + 1)
dmm_posterior <- extract_dmm_posterior(dmm_fit)
Z_dmm <- max.col(dmm_posterior, ties.method = "first")
t1 <- Sys.time()
elapsed_DMM <- as.numeric(difftime(t1, t0, units = "secs"))

mclust::adjustedRandIndex(Z_true, Z_dmm)

# run cor-tree ----------------------------------------------------------

# removed err_precision to prevent singularity
tree_depth <- set_tree_depth
cutoff_layer <- 4    # Layer for correlated nodes cutoff
warm_start <- 0
burnin <- 100       # Burn-in period
total_iter <- burnin + 50   # Total iterations
cov_interval <- 5
c_sigma2_vec <- 10  # Hyperparameter for sigma^2 vector. 1/c is the mean
sigma_mu2 <- 0.1     # Hyperparameter for sigma_mu^2
# Call the CorTree_sampler function
X_rowsum_raw = rowSums(X)
X_rowsum = dplyr::ntile(X_rowsum_raw, n_clus + 1)

cortree <- CorTree_sampler(X, 
                            # init_Z = Z_pam-1,
                            init_Z = X_rowsum-1,
                            n_clus = n_clus+1, 
                            tree_depth = tree_depth, 
                            cutoff_layer = cutoff_layer, 
                            total_iter = total_iter, 
                            burnin = burnin, 
                            cov_interval = cov_interval,
                            c_sigma2_vec = c_sigma2_vec, 
                            sigma_mu2,
                            warm_start = warm_start,
                            all_ind = F)
cortree$elapsed
# par(mfrow = c(1,1))
# plot(cortree$mcmc$loglik,type="l")

# clustering pattern
Z_cortree = apply(cortree$mcmc$Z[,-c(total_iter-burnin)],1,mean);table(Z_cortree)
Z_cortree = round(Z_cortree)

mclust::adjustedRandIndex(Z_true, Z_cortree)

# run ind-tree ----------------------------------------------------------

# removed err_precision to prevent singularity
tree_depth <- set_tree_depth
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
                            warm_start = warm_start,
                            all_ind = T)
indtree$elapsed
par(mfrow = c(1,1))
plot(indtree$mcmc$loglik,type="l")

# clustering pattern
Z_indtree = apply(indtree$mcmc$Z[,-c(total_iter-burnin)],1,mean);table(Z_indtree)
Z_indtree = round(Z_indtree); table(Z_indtree)

mclust::adjustedRandIndex(Z_true, Z_indtree)




# summarize results ----------------------------------------------------------

all_Z = cbind(Z_true, Z_kmeans, Z_pam, Z_dmm, Z_indtree, Z_cortree)
colnames(all_Z) = c("True","K-means","PAM","DMM","Ind-tree","Cor-tree")

all_ARI = cbind(
  mclust::adjustedRandIndex(Z_true, Z_kmeans),
  mclust::adjustedRandIndex(Z_true, Z_pam),
  mclust::adjustedRandIndex(Z_true, Z_dmm),
  mclust::adjustedRandIndex(Z_true, Z_indtree),
  mclust::adjustedRandIndex(Z_true, Z_cortree)
)
colnames(all_ARI) = colnames(all_Z)[-1]

# add time 
all_time = cbind(elapsed_knn, elapsed_PAM, elapsed_DMM, indtree$elapsed, cortree$elapsed)
colnames(all_time) = c("K-means","PAM","DMM","Ind-tree","Cor-tree")

out = list(
  all_Z = all_Z,
  all_ARI = all_ARI,
  all_time = all_time,
  ari = all_ARI,
  time = all_time,
  dmm = list(
    fit = dmm_fit,
    posterior = dmm_posterior,
    cluster = Z_dmm,
    elapsed = elapsed_DMM
  )
)

all_ARI

# draw a violin plot to show ARI of different methods under various settings
saveRDS(out,file.path(out_path,paste0("Sim1_clus_n",n_sample,"_rep",i_rep,".rds")))

# Run summary mode with:
# Rscript ./test/Sim1_clus_cortree.R summary
