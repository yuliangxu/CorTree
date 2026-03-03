CorTree_sampler_randominit <- function(
  X,
  n_clus,
  tree_depth,
  cutoff_layer,
  total_iter,
  burnin,
  n_start = 5L,
  init_Z_list = NULL,
  warm_start = 0L,
  c_sigma2_vec = 1.0,
  sigma_mu2 = 1.0,
  all_ind = FALSE,
  cov_interval = 1L,
  save_phi_trace = TRUE,
  save_cluster_cor_trace = TRUE,
  discard_collapsed_at_burnin = TRUE,
  verbose = TRUE
) {
  n_start <- as.integer(n_start)
  if (is.na(n_start) || n_start <= 0L) {
    stop("`n_start` must be a positive integer.")
  }
  if (!is.matrix(X)) {
    stop("`X` must be a matrix.")
  }
  n <- nrow(X)
  if (burnin >= total_iter) {
    stop("`burnin` must be smaller than `total_iter`.")
  }

  normalize_init <- function(z, n, n_clus) {
    z <- as.integer(z)
    if (length(z) != n) {
      stop("Each init_Z must have length equal to nrow(X).")
    }
    if (anyNA(z)) {
      stop("Each init_Z must not contain NA.")
    }
    if (any(z < 0L) || any(z > (n_clus - 1L))) {
      stop("Each init_Z must be in {0, ..., n_clus - 1}.")
    }
    z
  }

  provided_inits <- list()
  if (!is.null(init_Z_list)) {
    if (is.list(init_Z_list)) {
      provided_inits <- lapply(init_Z_list, normalize_init, n = n, n_clus = n_clus)
    } else if (is.matrix(init_Z_list)) {
      if (nrow(init_Z_list) != n) {
        stop("If `init_Z_list` is a matrix, it must have nrow(X) rows.")
      }
      provided_inits <- lapply(seq_len(ncol(init_Z_list)), function(j) {
        normalize_init(init_Z_list[, j], n = n, n_clus = n_clus)
      })
    } else {
      provided_inits <- list(normalize_init(init_Z_list, n = n, n_clus = n_clus))
    }
  }

  if (length(provided_inits) > n_start) {
    warning("More init_Z provided than `n_start`; only the first n_start are used.")
    provided_inits <- provided_inits[seq_len(n_start)]
  }

  n_random <- n_start - length(provided_inits)
  random_inits <- vector("list", n_random)
  for (i in seq_len(n_random)) {
    random_inits[[i]] <- sample.int(n_clus, size = n, replace = TRUE) - 1L
  }
  all_inits <- c(provided_inits, random_inits)

  fit_list <- vector("list", n_start)
  pi_trace_list <- vector("list", n_start)
  score_vec <- rep(-Inf, n_start)
  elapsed_vec <- rep(NA_real_, n_start)
  collapsed_at_burnin <- rep(FALSE, n_start)

  post_idx <- (burnin + 1L):total_iter
  for (s in seq_len(n_start)) {
    if (isTRUE(verbose)) {
      cat(sprintf("====== running %d-th random init ======\n", s))
    }
    fit_s <- CorTree_sampler(
      X = X,
      n_clus = n_clus,
      tree_depth = tree_depth,
      cutoff_layer = cutoff_layer,
      total_iter = total_iter,
      burnin = burnin,
      warm_start = warm_start,
      init_Z = all_inits[[s]],
      c_sigma2_vec = c_sigma2_vec,
      sigma_mu2 = sigma_mu2,
      all_ind = all_ind,
      cov_interval = cov_interval,
      save_phi_trace = save_phi_trace,
      save_cluster_cor_trace = save_cluster_cor_trace
    )

    fit_list[[s]] <- fit_s
    pi_trace_list[[s]] <- fit_s$mcmc$pi
    elapsed_vec[s] <- fit_s$elapsed

    z_chain_s <- fit_s$mcmc$Z
    if (!is.null(z_chain_s) && ncol(z_chain_s) >= 1L) {
      z_burnin_end <- as.integer(z_chain_s[, 1L])
      collapsed_s <- length(unique(z_burnin_end)) <= 1L
      collapsed_at_burnin[s] <- collapsed_s
      if (isTRUE(discard_collapsed_at_burnin) && collapsed_s) {
        if (isTRUE(verbose)) {
          cat(sprintf("====== discard %d-th init: collapsed to a single cluster at burn-in end ======\n", s))
        }
        score_vec[s] <- -Inf
        next
      }
    }

    loglik_s <- fit_s$mcmc$loglik
    score_vec[s] <- max(loglik_s[post_idx], na.rm = TRUE)
  }

  valid_idx <- which(is.finite(score_vec))
  best_idx <- if (length(valid_idx) > 0L) valid_idx[which.max(score_vec[valid_idx])] else NA_integer_

  list(
    best_fit = if (is.na(best_idx)) NULL else fit_list[[best_idx]],
    all_fit = fit_list,
    best_start_id = best_idx,
    best_init_Z = if (is.na(best_idx)) NULL else all_inits[[best_idx]],
    best_postburnin_loglik = if (is.na(best_idx)) NA_real_ else score_vec[best_idx],
    start_postburnin_loglik = score_vec,
    collapsed_at_burnin = collapsed_at_burnin,
    init_Z_all = all_inits,
    pi_mcmc_all = pi_trace_list,
    elapsed_all = elapsed_vec
  )
}
