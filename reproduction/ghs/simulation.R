#!/usr/bin/env Rscript
# Portable matched Sim1 reproduction. All generated files belong outside this repo.
options(stringsAsFactors = FALSE)
Sys.setenv(TZ = "UTC")

usage <- function() cat(paste(
  "Rscript reproduction/ghs/simulation.R --output DIR --task 1..300|all [--all-priors] [--library DIR]",
  "Rscript reproduction/ghs/simulation.R --output DIR --smoke [--all-priors] [--library DIR]",
  "Rscript reproduction/ghs/simulation.R --output DIR --summarize",
  "Smoke outputs are isolated from the 300 production tasks. No figures are generated.",
  sep = "\n"), "\n")

parse_args <- function(args) {
  out <- list(output = NULL, task = NULL, library = NULL, smoke = FALSE,
    all_priors = FALSE, summarize = FALSE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("--help", "-h")) { usage(); quit(save = "no") }
    if (key %in% c("--smoke", "--all-priors", "--summarize")) {
      out[[gsub("-", "_", substring(key, 3L))]] <- TRUE
    } else if (key %in% c("--output", "--task", "--library")) {
      i <- i + 1L
      if (i > length(args) || startsWith(args[i], "--")) stop("Missing value for ", key)
      out[[substring(key, 3L)]] <- args[i]
    } else stop("Unknown argument: ", key)
    i <- i + 1L
  }
  if (is.null(out$output)) stop("--output is required. Use --help for examples.")
  if (out$summarize) {
    if (out$smoke || !is.null(out$task)) stop("--summarize cannot be combined with --task or --smoke.")
  } else {
    if (is.null(out$task) && out$smoke) out$task <- "1"
    if (is.null(out$task)) stop("Specify --task 1..300|all, or --smoke.")
    if (out$task != "all" && !grepl("^[0-9]+$", out$task)) stop("Invalid task id.")
    ids <- if (out$task == "all") 1:300 else as.integer(out$task)
    if (anyNA(ids) || any(ids < 1L | ids > 300L)) stop("Task ids must be 1..300.")
    if (out$smoke && !identical(ids, 1L)) stop("--smoke uses only task 1, with N=20.")
    out$ids <- ids
  }
  out
}

# Resolve existing ancestors too, so a symlink cannot put outputs inside the repo.
absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(getwd(), path)
  suffix <- character()
  while (!dir.exists(path)) {
    parent <- dirname(path)
    if (identical(parent, path) || file.exists(path)) stop("Invalid output directory: ", path)
    suffix <- c(basename(path), suffix)
    path <- parent
  }
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  for (part in suffix) resolved <- if (part == "..") dirname(resolved) else
    if (part == ".") resolved else file.path(resolved, part)
  resolved
}

write_rds <- function(object, path) {
  # A unique sibling also permits different task processes to initialize the
  # same identical root protocol concurrently. tempfile() does not use R's RNG.
  temporary <- tempfile(pattern = paste0(basename(path), ".tmp-"), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Could not finalize ", path)
}
write_csv <- function(object, path) utils::write.csv(object, path, row.names = FALSE, na = "")
stamp <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
md5 <- function(paths) unname(tools::md5sum(paths))

prior_settings <- function(n, all_priors) {
  arm <- c("flat", "uniform", "soft_s0p1", "strong_soft_s0p1")
  a <- c(0, 0, 4, n - n / 3)
  # Preserve the archived Python manifest's left-to-right multiplication;
  # a * (s0^2) / 2 can differ by one ULP and change an entire replayed chain.
  rho <- c(0, 0, 2 * 0.1 * 0.1, (n - n / 3) * 0.1 * 0.1 / 2)
  upper <- c(Inf, 100, Inf, Inf)
  hierarchy <- rep(FALSE, 4L)
  if (all_priors) {
    arm <- c(arm, "exponential", "bounded_4N3", "bounded_2N3", "soft_s0p3", "soft_s1p0",
      "hier_s0p1", "hier_s0p3", "hier_s1p0")
    a <- c(a, 0, 4 * n / 3, 2 * n / 3, rep(4, 5L))
    rho <- c(rho, log(2), 0, 0, 2 * 0.3 * 0.3, 2, 2 * 0.1 * 0.1, 2 * 0.3 * 0.3, 2)
    upper <- c(upper, Inf, 100, 100, rep(Inf, 5L))
    hierarchy <- c(hierarchy, rep(FALSE, 5L), rep(TRUE, 3L))
  }
  data.frame(arm = arm, a = a, rho = rho, upper = upper, hierarchy = hierarchy,
    c = 3, r = 4, s = 2,
    dahl = ifelse(grepl("^(soft_|hier_|strong_soft_)", arm), "integer_contingency", "package"))
}

# Historical soft-prior studies used this exact integer objective. The earlier
# studies used the package's matrix objective; retain each historical tie rule.
dahl_integer <- function(Z) {
  draws <- ncol(Z)
  codes <- lapply(seq_len(draws), function(i) match(Z[, i], unique(Z[, i])))
  K <- vapply(codes, max, integer(1L))
  overlap <- matrix(0, draws, draws)
  for (i in seq_len(draws)) {
    overlap[i, i] <- sum(as.double(tabulate(codes[[i]]))^2)
    if (i < draws) for (j in seq.int(i + 1L, draws)) {
      cells <- tabulate(codes[[i]] + K[i] * (codes[[j]] - 1L), K[i] * K[j])
      overlap[i, j] <- overlap[j, i] <- sum(as.double(cells)^2)
    }
  }
  objective <- draws * diag(overlap) - 2 * rowSums(overlap)
  best <- which.min(objective)
  list(Z = as.integer(Z[, best]), best_iter = best, loss = objective / draws + mean(overlap))
}

generate_input <- function(task_id, smoke) {
  paper_n <- c(200L, 400L, 600L)[(task_id - 1L) %/% 100L + 1L]
  n <- if (smoke) 20L else paper_n
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(2025L + task_id)
  labels_true <- sample(1:2, n, prob = c(0.6, 0.4), replace = TRUE)
  X <- matrix(0, nrow = n, ncol = 1000L)
  gen_X <- function(n, cluster_id) {
    weight <- stats::rbeta(1L, 10, 10)
    n1 <- floor(n * weight)
    n2 <- n - n1
    switch(cluster_id,
      c(stats::rbeta(n1, 2, 6), stats::rbeta(n2, 6, 2)),
      c(stats::rbeta(n1, 1, 1), stats::rbeta(n2, 3, 3)))
  }
  histogram_breaks <- seq(0, 1, length.out = 1001L)
  # Keep this group loop and random-call ordering identical to the original.
  for (cluster_id in 1:2) {
    idx <- which(labels_true == cluster_id)
    if (!length(idx)) next
    samples <- replicate(length(idx), gen_X(sample(1000:5000, 1L), cluster_id), simplify = FALSE)
    X[idx, ] <- do.call(rbind, lapply(samples, function(x)
      graphics::hist(x, breaks = histogram_breaks, plot = FALSE)$counts))
  }
  X <- matrix(as.integer(round(X)), nrow = nrow(X), ncol = ncol(X),
    dimnames = list(paste0("sample_", seq_len(nrow(X))), paste0("taxon_", seq_len(ncol(X)))))
  initial_labels <- dplyr::ntile(rowSums(X), 3L) - 1L
  list(X = X, labels_true = labels_true, initial_labels = initial_labels,
    RNGkind = RNGkind(), rng_after_data = .Random.seed, paper_n = paper_n)
}

fit_diagnostics <- function(fit, spec, n, total_iter, burnin, initial_labels, independent = FALSE) {
  Z <- as.matrix(fit$mcmc$Z)
  retained <- total_iter - burnin
  stopifnot(identical(dim(Z), c(as.integer(n), as.integer(retained))),
    all(is.finite(Z)), all(Z %in% 0:2))
  weights <- as.matrix(fit$mcmc$pi)
  stopifnot(identical(dim(weights), c(3L, as.integer(total_iter))),
    all(is.finite(weights)), all(weights >= 0), max(abs(colSums(weights) - 1)) < 1e-8,
    length(fit$mcmc$loglik) == total_iter, all(is.finite(fit$mcmc$loglik)))
  occupancy <- vapply(seq_len(retained), function(i) tabulate(Z[, i] + 1L, 3L), integer(3L))
  canonical <- vapply(seq_len(retained), function(i) match(Z[, i], unique(Z[, i])), integer(n))
  adjacent_ari <- vapply(seq.int(2L, retained), function(i)
    mclust::adjustedRandIndex(Z[, i - 1L], Z[, i]), numeric(1L))
  same <- colSums(canonical[, -1L, drop = FALSE] != canonical[, -retained, drop = FALSE]) == 0L
  mixing <- data.frame(draw = seq_len(retained), iteration = burnin + seq_len(retained),
    occupied_clusters = colSums(occupancy > 0), adjacent_partition_ari = c(NA, adjacent_ari),
    adjacent_identical = c(NA, same), raw_label_change_fraction = c(NA,
      colMeans(Z[, -1L, drop = FALSE] != Z[, -retained, drop = FALSE])),
    initial_partition_ari = vapply(seq_len(retained), function(i)
      mclust::adjustedRandIndex(initial_labels, Z[, i]), numeric(1L)))
  guard <- fit$mcmc$covariance_safeguards
  stopifnot(is.list(guard), all(c("enabled", "regularizations", "updates_skipped") %in% names(guard)))
  precision <- data.frame(draw = integer(), iteration = integer(), component = integer(),
    occupancy = integer(), min_eigenvalue = numeric(), max_diagonal = numeric(),
    condition_number = numeric(), log_determinant = numeric(), scale_t = numeric(), shared_b = numeric())
  if (!independent) {
    prior <- fit$mcmc$ghs_prior
    proper <- spec$arm != "flat"
    stopifnot(isTRUE(prior$active), identical(isTRUE(prior$proper), proper),
      isTRUE(all.equal(prior$diag_rate, spec$rho)), isTRUE(all.equal(prior$diag_upper, spec$upper)),
      isTRUE(all.equal(prior$det_df, spec$a)), length(fit$mcmc$Sigma_inv) == retained)
    if (proper) stopifnot(!isTRUE(guard$enabled), sum(guard$regularizations) == 0,
      sum(guard$updates_skipped) == 0)
    hierarchy <- fit$mcmc$ghs_scale
    if (spec$hierarchy) stopifnot(isTRUE(hierarchy$active),
      identical(dim(as.matrix(hierarchy$t)), c(3L, as.integer(retained))),
      length(hierarchy$b) == retained, all(is.finite(hierarchy$t)), all(hierarchy$t > 0),
      all(is.finite(hierarchy$b)), all(hierarchy$b > 0), hierarchy$shape == spec$c,
      hierarchy$rate_shape == spec$r, hierarchy$rate_rate == spec$s)
    else if (!is.null(hierarchy)) stopifnot(!isTRUE(hierarchy$active))
    precision <- do.call(rbind, lapply(seq_len(retained), function(i) {
      cube <- fit$mcmc$Sigma_inv[[i]]
      stopifnot(identical(dim(cube), c(31L, 31L, 3L)), all(is.finite(cube)))
      do.call(rbind, lapply(1:3, function(k) {
        omega <- cube[, , k]
        ev <- eigen((omega + t(omega)) / 2, symmetric = TRUE, only.values = TRUE)$values
        stopifnot(max(abs(omega - t(omega))) < 1e-8 * max(1, max(abs(omega))),
          all(is.finite(ev)), min(ev) > 0, all(diag(omega) > 0))
        if (is.finite(spec$upper)) stopifnot(all(diag(omega) < spec$upper))
        data.frame(draw = i, iteration = burnin + i, component = k - 1L, occupancy = occupancy[k, i],
          min_eigenvalue = min(ev), max_diagonal = max(diag(omega)),
          condition_number = max(ev) / min(ev), log_determinant = sum(log(ev)),
          scale_t = if (spec$hierarchy) hierarchy$t[k, i] else 1,
          shared_b = if (spec$hierarchy) hierarchy$b[i] else NA_real_)
      }))
    }))
  }
  list(mixing = mixing, precision = precision, occupancy = occupancy,
    summary = list(retained_unique_partitions = nrow(unique(t(canonical))),
      fraction_adjacent_identical = mean(same), mean_adjacent_partition_ari = mean(adjacent_ari),
      min_occupied_clusters = min(mixing$occupied_clusters), max_occupied_clusters = max(mixing$occupied_clusters),
      min_precision_eigenvalue = if (nrow(precision)) min(precision$min_eigenvalue) else NA_real_,
      max_precision_diagonal = if (nrow(precision)) max(precision$max_diagonal) else NA_real_,
      regularizations = sum(guard$regularizations), updates_skipped = sum(guard$updates_skipped),
      empty_precision_updates = sum(fit$mcmc$empty_precision_updates), finite_loglik = TRUE))
}

validate_complete <- function(path, expected_methods = NULL) {
  marker <- file.path(path, "COMPLETE")
  manifest <- file.path(path, "completion.rds")
  if (!file.exists(marker) || !file.exists(manifest)) stop("Missing completion marker: ", path)
  stopifnot(identical(readLines(marker, warn = FALSE), md5(manifest)))
  complete <- readRDS(manifest)
  stopifnot(identical(md5(file.path(path, names(complete$artifact_md5))), unname(complete$artifact_md5)))
  metrics <- utils::read.csv(file.path(path, "metrics.csv"), stringsAsFactors = FALSE)
  stopifnot(all(is.finite(metrics$ari)), !anyDuplicated(metrics$method),
    identical(sort(metrics$method), sort(complete$methods)))
  if (!is.null(expected_methods)) stopifnot(identical(sort(metrics$method), sort(expected_methods)))
  list(completion = complete, metrics = metrics)
}

run_task <- function(task_id, opt, protocol, root) {
  out <- if (opt$smoke) file.path(root, "smoke", "task_001") else
    file.path(root, "simulation", sprintf("task_%03d", task_id))
  if (file.exists(file.path(out, "COMPLETE"))) {
    previous <- validate_complete(out, protocol$methods)
    stopifnot(identical(previous$completion$protocol, protocol))
    message("Already verified complete: ", out)
    return(invisible(NULL))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(out, ".running")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Task is locked: ", lock,
    ". Remove this directory only after confirming that its process has stopped.")
  writeLines(c(paste("pid", Sys.getpid()), stamp()), file.path(lock, "owner.txt"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  old_wd <- getwd()
  setwd(out) # Some baseline packages may write diagnostic files in the cwd.
  on.exit(setwd(old_wd), add = TRUE)
  stage <- "generate_input"
  tryCatch({
    input <- generate_input(task_id, opt$smoke)
    n <- nrow(input$X)
    total_iter <- if (opt$smoke) 6L else 150L
    burnin <- if (opt$smoke) 2L else 100L
    specs <- prior_settings(n, opt$all_priors)
    metadata <- list(task_id = task_id, replicate_id = (task_id - 1L) %% 100L + 1L,
      seed = 2025L + task_id, n = n, paper_n = input$paper_n, smoke = opt$smoke,
      total_iter = total_iter, burnin = burnin, fitted_K = 3L, tree_depth = 6L,
      cortree_cutoff = 4L, indtree_cutoff = 3L, cov_interval = 5L, c_sigma2_vec = 10,
      sigma_mu2 = 0.1, warm_start = 0L, prior_settings = specs, protocol = protocol,
      started_utc = stamp(), package_library = find.package("CorTree"),
      rng_protocol = "Kmeans -> PAM -> DMM -> flat CorTree -> IndTree; reset to flat entry for each proper arm",
      generator_note = "Group loop with floor(m*Beta(10,10)) first-component samples; this is the historical code generator.")
    input$metadata <- metadata
    write_rds(input, "input.rds")
    write_rds(metadata, "metadata.rds")
    writeLines(capture.output(utils::sessionInfo()), "sessionInfo.txt")
    write_csv(specs, "prior_settings.csv")
    metrics <- list()
    capture_warnings <- function(fun) {
      messages <- character()
      start <- proc.time()[["elapsed"]]
      value <- withCallingHandlers(fun(), warning = function(w) {
        messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
      })
      list(value = value, warnings = messages, wall_sec = proc.time()[["elapsed"]] - start)
    }
    add_metric <- function(method, labels, elapsed, warnings, diagnostic = NULL, spec = NULL) {
      row <- data.frame(task_id = task_id, replicate_id = metadata$replicate_id, smoke = opt$smoke,
        n = n, seed = metadata$seed, method = method, ari = mclust::adjustedRandIndex(input$labels_true, labels),
        estimated_clusters = length(unique(labels)), elapsed_sec = elapsed,
        warnings = paste(warnings, collapse = " | "), a = if (is.null(spec)) NA_real_ else spec$a,
        rho = if (is.null(spec)) NA_real_ else spec$rho,
        upper = if (is.null(spec)) NA_real_ else spec$upper,
        hierarchy = if (is.null(spec)) FALSE else spec$hierarchy,
        retained_unique_partitions = if (is.null(diagnostic)) NA_integer_ else diagnostic$retained_unique_partitions,
        fraction_adjacent_identical = if (is.null(diagnostic)) NA_real_ else diagnostic$fraction_adjacent_identical,
        regularizations = if (is.null(diagnostic)) 0 else diagnostic$regularizations,
        updates_skipped = if (is.null(diagnostic)) 0 else diagnostic$updates_skipped,
        min_precision_eigenvalue = if (is.null(diagnostic)) NA_real_ else diagnostic$min_precision_eigenvalue,
        max_precision_diagonal = if (is.null(diagnostic)) NA_real_ else diagnostic$max_precision_diagonal)
      stopifnot(is.finite(row$ari), is.finite(row$elapsed_sec))
      metrics[[method]] <<- row
      write_csv(do.call(rbind, metrics), "metrics.partial.csv")
    }
    # These calls, including dmn's DEFAULT seed, must precede the saved RNG.
    stage <- "Kmeans"
    km <- capture_warnings(function() stats::kmeans(scale(input$X), centers = 3L, nstart = 25L))
    write_rds(km, "Kmeans.rds")
    add_metric("Kmeans", km$value$cluster, km$wall_sec, km$warnings)
    stage <- "PAM"
    pam <- capture_warnings(function() cluster::pam(input$X, 3L))
    write_rds(pam, "PAM.rds")
    add_metric("PAM", pam$value$clustering, pam$wall_sec, pam$warnings)
    stage <- "DMM"
    dmm <- capture_warnings(function() DirichletMultinomial::dmn(input$X, k = 3L))
    posterior <- as.matrix(DirichletMultinomial::mixture(dmm$value))
    if (nrow(posterior) != n && ncol(posterior) == n) posterior <- t(posterior)
    stopifnot(nrow(posterior) == n, ncol(posterior) == 3L, all(is.finite(posterior)))
    write_rds(list(fit = dmm$value, posterior = posterior, warnings = dmm$warnings), "DMM.rds")
    add_metric("DMM", max.col(posterior, ties.method = "first"), dmm$wall_sec, dmm$warnings)
    rng_cortree <- .Random.seed
    write_rds(rng_cortree, "rng_before_cortree.rds")
    run_sampler <- function(spec, independent = FALSE) {
      method <- if (independent) "IndTree" else spec$arm
      stage <<- method
      message("task ", task_id, " / ", method, " / N=", n, if (opt$smoke) " (SMOKE)" else "")
      arm_dir <- file.path(out, method)
      dir.create(arm_dir, showWarnings = FALSE)
      write_rds(.Random.seed, file.path(arm_dir, "rng_before_sampler.rds"))
      args <- list(X = input$X, n_clus = 3L, tree_depth = 6L,
        cutoff_layer = if (independent) 3L else 4L, total_iter = total_iter, burnin = burnin,
        warm_start = 0L, init_Z = input$initial_labels, c_sigma2_vec = 10, sigma_mu2 = 0.1,
        all_ind = independent, cov_interval = 5L, save_phi_trace = FALSE, save_cluster_cor_trace = FALSE)
      # Leave the original flat/IndTree calls untouched, including default args.
      if (!independent && spec$arm != "flat") args <- c(args, list(ghs_diag_rate = spec$rho,
        ghs_diag_upper = spec$upper, ghs_det_df = spec$a, ghs_scale_hierarchy = spec$hierarchy,
        ghs_scale_shape = spec$c, ghs_scale_rate_shape = spec$r, ghs_scale_rate_rate = spec$s))
      result <- capture_warnings(function() do.call(CorTree::CorTree_sampler, args))
      write_rds(.Random.seed, file.path(arm_dir, "rng_after_sampler.rds"))
      fit <- result$value
      diagnostic <- fit_diagnostics(fit, spec, n, total_iter, burnin, input$initial_labels, independent)
      dahl <- if (!independent && spec$dahl == "integer_contingency") dahl_integer(as.matrix(fit$mcmc$Z)) else
        CorTree::dahl_clustering(as.matrix(fit$mcmc$Z))
      fit$dahl <- dahl
      fit$metadata <- metadata
      fit$prior_setting <- if (independent) NULL else spec
      fit$warnings <- result$warnings
      write_rds(fit, file.path(arm_dir, "fit.rds"))
      write_rds(diagnostic, file.path(arm_dir, "diagnostics.rds"))
      write_csv(diagnostic$mixing, file.path(arm_dir, "mixing.csv"))
      write_csv(diagnostic$precision, file.path(arm_dir, "precision.csv"))
      write_csv(data.frame(sample_id = rownames(input$X), truth = input$labels_true,
        initial = input$initial_labels, dahl = dahl$Z), file.path(arm_dir, "assignments.csv"))
      writeLines(result$warnings, file.path(arm_dir, "warnings.txt"))
      if (!independent && spec$arm != "flat" && length(result$warnings))
        stop("Proper-prior sampler emitted warnings; retained files require inspection.")
      add_metric(method, dahl$Z, if (is.null(fit$elapsed)) result$wall_sec else as.numeric(fit$elapsed),
        result$warnings, diagnostic$summary, if (independent) NULL else spec)
      invisible(NULL)
    }
    run_sampler(specs[1L, ])
    # Preserve the old corrected workflow's IndTree RNG, immediately after flat.
    write_rds(.Random.seed, "rng_before_indtree.rds")
    run_sampler(specs[1L, ], independent = TRUE)
    for (i in seq.int(2L, nrow(specs))) {
      assign(".Random.seed", rng_cortree, envir = .GlobalEnv)
      run_sampler(specs[i, ])
    }
    table <- do.call(rbind, metrics)
    stopifnot(identical(sort(table$method), sort(protocol$methods)))
    write_csv(table, "metrics.csv")
    # Hash only this workflow's artifacts; COMPLETE is written last.
    files <- list.files(out, recursive = TRUE, all.files = FALSE, full.names = FALSE)
    files <- files[!files %in% c("completion.rds", "COMPLETE", "FAILED.txt", "metrics.partial.csv") & !grepl("\\.tmp$", files)]
    completion <- list(task_id = task_id, smoke = opt$smoke, methods = protocol$methods,
      completed_utc = stamp(), protocol = protocol, artifact_md5 = setNames(md5(file.path(out, files)), files))
    stopifnot(!anyNA(completion$artifact_md5))
    write_rds(completion, "completion.rds")
    writeLines(md5("completion.rds"), "COMPLETE")
    if (file.exists("FAILED.txt")) unlink("FAILED.txt")
    validate_complete(out, protocol$methods)
    message("Completed: ", out)
  }, error = function(e) {
    writeLines(c(stamp(), paste("stage:", stage), conditionMessage(e)), file.path(out, "FAILED.txt"))
    stop(e)
  })
  invisible(NULL)
}

summarize_results <- function(root) {
  protocol_path <- file.path(root, "protocol.rds")
  if (!file.exists(protocol_path)) stop("No production protocol.rds found; smoke runs cannot be summarized as production.")
  protocol <- readRDS(protocol_path)
  metrics <- list()
  status <- data.frame(task_id = 1:300, n = rep(c(200L, 400L, 600L), each = 100L),
    state = "missing", detail = "")
  for (task in 1:300) {
    path <- file.path(root, "simulation", sprintf("task_%03d", task))
    if (!file.exists(file.path(path, "COMPLETE"))) {
      if (file.exists(file.path(path, "FAILED.txt"))) status$state[task] <- "failed"
      else if (dir.exists(path)) status$state[task] <- "incomplete"
      next
    }
    item <- tryCatch({
      item <- validate_complete(path, protocol$methods)
      stopifnot(identical(item$completion$protocol, protocol), item$completion$task_id == task,
        !item$completion$smoke, all(item$metrics$task_id == task), !any(item$metrics$smoke),
        all(item$metrics$n == status$n[task]), all(item$metrics$seed == 2025L + task))
      item
    }, error = function(e) e)
    if (inherits(item, "error")) {
      status$state[task] <- "invalid"
      status$detail[task] <- conditionMessage(item)
    } else {
      status$state[task] <- "complete"
      metrics[[as.character(task)]] <- item$metrics
    }
  }
  all_metrics <- if (length(metrics)) do.call(rbind, metrics) else data.frame()
  summary <- do.call(rbind, lapply(c(200L, 400L, 600L), function(n) {
    do.call(rbind, lapply(protocol$methods, function(method) {
      rows <- if (nrow(all_metrics)) all_metrics[all_metrics$n == n & all_metrics$method == method, ] else data.frame()
      count <- nrow(rows)
      complete <- count == 100L
      differences <- numeric()
      if (complete) {
        flat <- all_metrics[all_metrics$n == n & all_metrics$method == "flat", ]
        differences <- rows$ari - flat$ari[match(rows$task_id, flat$task_id)]
        stopifnot(length(differences) == 100L, all(is.finite(differences)))
      }
      data.frame(n = n, method = method, complete_replicates = count, expected_replicates = 100L,
        status = if (complete) "complete" else "INCOMPLETE",
        mean_ari = if (complete) mean(rows$ari) else NA_real_,
        sd_ari = if (complete) stats::sd(rows$ari) else NA_real_,
        partial_mean_ari = if (count && !complete) mean(rows$ari) else NA_real_,
        partial_sd_ari = if (count > 1L && !complete) stats::sd(rows$ari) else NA_real_,
        mean_paired_difference_vs_flat = if (complete) mean(differences) else NA_real_,
        paired_mcse = if (complete) stats::sd(differences) / sqrt(100) else NA_real_,
        invariant_retained_fits = if (count) sum(rows$retained_unique_partitions == 1L, na.rm = TRUE) else 0L,
        fits_with_warnings = if (count) sum(!is.na(rows$warnings) & nzchar(rows$warnings)) else 0L)
    }))
  }))
  destination <- file.path(root, "summary")
  dir.create(destination, showWarnings = FALSE)
  write_csv(status, file.path(destination, "task_status.csv"))
  write_csv(all_metrics, file.path(destination, "metrics.csv"))
  write_csv(summary, file.path(destination, "simulation_summary.csv"))
  final <- all(status$state == "complete")
  writeLines(c(if (final) "COMPLETE: 300/300 matched simulation datasets." else
    sprintf("INCOMPLETE: %d/300 matched simulation datasets validated.", sum(status$state == "complete")),
    "Each final N/method mean requires all 100 matched replicates; partial means have separate columns.",
    "ARI uses Dahl partitions; an invariant retained chain is a mixing warning, not proof of convergence.",
    "Paired MCSE is sd(method ARI - matched flat ARI)/sqrt(100), not a posterior uncertainty estimate.",
    "These are corrected reruns, not the original published table or independent tuning holdouts."),
    file.path(destination, "STATUS.txt"))
  print(summary, row.names = FALSE)
  message(readLines(file.path(destination, "STATUS.txt"), n = 1L))
  invisible(summary)
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) != 1L) stop("Run this file with Rscript.")
  script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  repo <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
  root <- absolute_path(opt$output)
  if (identical(root, repo) || startsWith(root, paste0(repo, "/"))) stop("--output must be outside the source repository.")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  if (opt$summarize) { summarize_results(root); return(invisible(NULL)) }
  if (!is.null(opt$library)) .libPaths(c(normalizePath(opt$library, mustWork = TRUE), .libPaths()))
  required <- c("CorTree", "dplyr", "cluster", "mclust", "DirichletMultinomial")
  for (package in required) if (!requireNamespace(package, quietly = TRUE))
    stop("Install required package: ", package, ". See SIMULATION.md.")
  if (!is.null(opt$library)) stopifnot(identical(normalizePath(find.package("CorTree")),
    normalizePath(file.path(opt$library, "CorTree"), mustWork = TRUE)))
  stopifnot(all(c("ghs_det_df", "ghs_diag_upper", "ghs_scale_hierarchy", "ghs_scale_rate_shape") %in%
    names(formals(CorTree::CorTree_sampler))))
  library_files <- list.files(find.package("CorTree"), pattern = "\\.(so|dll|dylib)$", recursive = TRUE, full.names = TRUE)
  stopifnot(length(library_files) > 0L, !anyNA(md5(library_files)))
  methods <- c("Kmeans", "PAM", "DMM", "flat", "IndTree", prior_settings(200, opt$all_priors)$arm[-1L])
  protocol <- list(schema_version = 1L, all_priors = opt$all_priors, methods = methods,
    smoke = opt$smoke, R_version = R.version.string,
    package_versions = setNames(vapply(required, function(p) as.character(utils::packageVersion(p)), character(1L)), required),
    sampler_binary_md5 = setNames(md5(library_files), basename(library_files)), runner_md5 = md5(script),
    RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"), seed_rule = "2025 + task_id",
    production_tasks = 1:300, production_n = rep(c(200L, 400L, 600L), each = 100L))
  location <- if (opt$smoke) file.path(root, "smoke") else root
  dir.create(location, recursive = TRUE, showWarnings = FALSE)
  protocol_path <- file.path(location, "protocol.rds")
  if (file.exists(protocol_path)) {
    if (!identical(readRDS(protocol_path), protocol)) stop("Output protocol differs (settings, runner, R, or package). Use a fresh output directory.")
  } else write_rds(protocol, protocol_path)
  for (task_id in opt$ids) run_task(task_id, opt, protocol, root)
  # Individual array tasks must not race when writing shared summary files.
  if (!opt$smoke && identical(opt$task, "all")) summarize_results(root)
}

main()
