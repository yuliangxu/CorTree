CorTree_sampler_randominit_heldout <- function(
  X,
  train_idx,
  test_idx,
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
  save_phi_trace = FALSE,
  save_cluster_cor_trace = FALSE,
  n_phi_mc = 1L,
  discard_collapsed_at_burnin = TRUE,
  refit_full_data = TRUE,
  verbose = TRUE
) {
  n_start <- as.integer(n_start)
  if (is.na(n_start) || n_start <= 0L) {
    stop("`n_start` must be a positive integer.")
  }
  if (anyDuplicated(c(train_idx, test_idx)) > 0L) {
    stop("`train_idx` and `test_idx` must be disjoint.")
  }
  if (length(train_idx) == 0L || length(test_idx) == 0L) {
    stop("`train_idx` and `test_idx` must be non-empty.")
  }

  X_train <- X[train_idx, , drop = FALSE]
  X_test <- X[test_idx, , drop = FALSE]

  train_n <- nrow(X_train)
  n <- nrow(X)
  normalize_init <- function(z, n_expected, n_clus) {
    z <- as.integer(z)
    if (length(z) != n_expected) {
      stop("Each init_Z must have length equal to nrow(X) or nrow(X_train).")
    }
    if (anyNA(z)) {
      stop("Each init_Z must not contain NA.")
    }
    if (any(z < 0L) || any(z > (n_clus - 1L))) {
      stop("Each init_Z must be in {0, ..., n_clus - 1}.")
    }
    z
  }
  to_full_init <- function(z) {
    z <- as.integer(z)
    if (length(z) == n) {
      return(normalize_init(z, n_expected = n, n_clus = n_clus))
    }
    if (length(z) == train_n) {
      z_train <- normalize_init(z, n_expected = train_n, n_clus = n_clus)
      z_full <- sample.int(n_clus, size = n, replace = TRUE) - 1L
      z_full[train_idx] <- z_train
      return(z_full)
    }
    stop("Each provided init_Z must have length nrow(X) or nrow(X_train).")
  }

  provided_full_inits <- list()
  provided_names <- character(0)
  if (!is.null(init_Z_list)) {
    if (is.list(init_Z_list)) {
      provided_full_inits <- lapply(init_Z_list, to_full_init)
      provided_names <- names(init_Z_list)
    } else {
      provided_full_inits <- list(to_full_init(init_Z_list))
      provided_names <- names(provided_full_inits)
    }
  }

  if (length(provided_full_inits) > n_start) {
    warning("More init_Z provided than `n_start`; only the first n_start are used.")
    provided_full_inits <- provided_full_inits[seq_len(n_start)]
    if (length(provided_names) > 0L) {
      provided_names <- provided_names[seq_len(n_start)]
    }
  }

  n_random <- n_start - length(provided_full_inits)
  random_full_inits <- vector("list", n_random)
  for (i in seq_len(n_random)) {
    random_full_inits[[i]] <- sample.int(n_clus, size = n, replace = TRUE) - 1L
  }
  all_init_Z_full <- c(provided_full_inits, random_full_inits)
  provided_names <- if (length(provided_names) == length(provided_full_inits)) provided_names else rep("", length(provided_full_inits))
  random_names <- if (n_random > 0L) paste0("random_", seq_len(n_random)) else character(0)
  all_names <- c(provided_names, random_names)
  if (length(all_names) == length(all_init_Z_full)) {
    names(all_init_Z_full) <- all_names
  }

  init_Z_list_train <- lapply(all_init_Z_full, function(z_full) z_full[train_idx])
  names(init_Z_list_train) <- names(all_init_Z_full)

  fit_multi <- CorTree_sampler_randominit(
    X = X_train,
    n_clus = n_clus,
    tree_depth = tree_depth,
    cutoff_layer = cutoff_layer,
    total_iter = total_iter,
    burnin = burnin,
    n_start = n_start,
    init_Z_list = init_Z_list_train,
    warm_start = warm_start,
    c_sigma2_vec = c_sigma2_vec,
    sigma_mu2 = sigma_mu2,
    all_ind = all_ind,
    cov_interval = cov_interval,
    save_phi_trace = save_phi_trace,
    save_cluster_cor_trace = save_cluster_cor_trace,
    discard_collapsed_at_burnin = discard_collapsed_at_burnin,
    verbose = verbose
  )

  hlpd_by_start <- rep(NA_real_, length(fit_multi$all_fit))
  hlpd_elapsed_by_start <- rep(NA_real_, length(fit_multi$all_fit))
  pred_detail <- vector("list", length(fit_multi$all_fit))
  for (s in seq_along(fit_multi$all_fit)) {
    if (isTRUE(discard_collapsed_at_burnin) && isTRUE(fit_multi$collapsed_at_burnin[s])) {
      if (isTRUE(verbose)) {
        cat(sprintf("====== skip HLPD for %d-th random init: discarded due to burn-in collapse ======\n", s))
      }
      next
    }
    if (isTRUE(verbose)) {
      cat(sprintf("====== evaluating HLPD for %d-th random init ======\n", s))
    }
    t_hlpd <- system.time({
      pred_s <- CorTree_heldout_logpred(
        X_test = X_test,
        mcmc = fit_multi$all_fit[[s]]$mcmc,
        tree_depth = tree_depth,
        cutoff_layer = cutoff_layer,
        burnin = burnin,
        all_ind = all_ind,
        n_phi_mc = n_phi_mc
      )
    })
    hlpd_elapsed_by_start[s] <- unname(t_hlpd[["elapsed"]])
    if (isTRUE(verbose)) {
      cat(sprintf("====== HLPD eval time for %d-th random init: %.3f sec ======\n", s, hlpd_elapsed_by_start[s]))
    }
    hlpd_by_start[s] <- pred_s$hlpd
    pred_detail[[s]] <- pred_s
  }

  valid_idx <- which(is.finite(hlpd_by_start))
  if (length(valid_idx) == 0L) {
    stop("All starts were discarded or had invalid HLPD; no valid held-out fit to select.")
  }
  best_idx <- valid_idx[which.max(hlpd_by_start[valid_idx])]
  best_fit_heldout <- fit_multi$all_fit[[best_idx]]
  best_init_Z_heldout <- fit_multi$init_Z_all[[best_idx]]
  best_init_Z_full <- all_init_Z_full[[best_idx]]

  full_fit <- NULL
  full_fit_elapsed <- NA_real_
  if (isTRUE(refit_full_data)) {
    if (isTRUE(verbose)) {
      cat(sprintf("====== refitting on full data using selected %d-th init_Z ======\n", best_idx))
    }
    t_full <- system.time({
      full_fit <- CorTree_sampler(
        X = X,
        n_clus = n_clus,
        tree_depth = tree_depth,
        cutoff_layer = cutoff_layer,
        total_iter = total_iter,
        burnin = burnin,
        warm_start = warm_start,
        init_Z = best_init_Z_full,
        c_sigma2_vec = c_sigma2_vec,
        sigma_mu2 = sigma_mu2,
        all_ind = all_ind,
        cov_interval = cov_interval,
        save_phi_trace = save_phi_trace,
        save_cluster_cor_trace = save_cluster_cor_trace
      )
    })
    full_fit_elapsed <- unname(t_full[["elapsed"]])
    if (isTRUE(verbose)) {
      cat(sprintf("====== full-data refit time: %.3f sec ======\n", full_fit_elapsed))
    }
  }

  list(
    best_start_id_heldout = best_idx,
    best_fit_heldout = best_fit_heldout,
    best_init_Z_heldout = best_init_Z_heldout,
    best_init_Z_full = best_init_Z_full,
    heldout_hlpd_by_start = hlpd_by_start,
    train_elapsed_by_start = fit_multi$elapsed_all,
    discarded_collapsed_by_start = fit_multi$collapsed_at_burnin,
    heldout_hlpd_elapsed_by_start = hlpd_elapsed_by_start,
    heldout_pred_detail = pred_detail,
    init_Z_full_all = all_init_Z_full,
    full_fit_from_combined_init = full_fit,
    full_fit_from_best_start_init = full_fit,
    full_fit_elapsed = full_fit_elapsed,
    fit_multi = fit_multi,
    train_idx = train_idx,
    test_idx = test_idx
  )
}
