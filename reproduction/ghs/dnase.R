#!/usr/bin/env Rscript
# Standalone reproduction; input data and generated artifacts stay outside git.
usage <- function() cat(paste0(
  "Usage: Rscript reproduction/ghs/dnase.R --output DIR [options]\n",
  "  --data-dir DIR       Input RDS directory, or set CORTREE_DNA_DATA\n",
  "  --dataset REST|NRF1  Default REST; one dataset per command\n",
  "  --init all|rowsum|pwm|tss_dist   Default all\n",
  "  --prior selected|all|NAME       Default selected; see DNASE.md\n",
  "  --all-priors         Alias for --prior all\n",
  "  --library DIR        Optional installed R library containing CorTree\n",
  "  --smoke              First 64 rows, 6 iterations/2 burn-in; not report results\n",
  "  --help               Show help without fitting\n"))

script_path <- function() {
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argument) != 1L) stop("Run this script with Rscript.")
  normalizePath(sub("^--file=", "", argument), mustWork = TRUE)
}

# Resolve an output that does not yet exist, including symlinked ancestors.
absolute_output <- function(path) {
  candidate <- path.expand(path)
  if (!grepl("^/|^[A-Za-z]:[/\\\\]", candidate)) candidate <- file.path(getwd(), candidate)
  suffix <- character()
  while (!dir.exists(candidate)) {
    parent <- dirname(candidate)
    if (identical(parent, candidate)) stop("Cannot resolve output directory.")
    suffix <- c(basename(candidate), suffix)
    candidate <- parent
  }
  candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  for (part in suffix) {
    if (part == "..") candidate <- dirname(candidate)
    else if (part != ".") candidate <- file.path(candidate, part)
  }
  candidate
}

read_options <- function() {
  args <- commandArgs(TRUE)
  if (!length(args) || "--help" %in% args) { usage(); quit(save = "no", status = 0L) }
  opts <- list(out = NULL, data_dir = Sys.getenv("CORTREE_DNA_DATA"),
    dataset = "REST", init = "all", prior = "selected", library = NULL, smoke = FALSE)
  i <- 1L
  while (i <= length(args)) {
    name <- sub("^--", "", args[i]); name <- gsub("-", "_", name, fixed = TRUE)
    if (name == "all_priors") { opts$prior <- "all"; i <- i + 1L; next }
    if (name == "output") name <- "out"
    if (!startsWith(args[i], "--") || !name %in% names(opts)) stop("Unknown option: ", args[i])
    if (name == "smoke") opts$smoke <- TRUE else {
      i <- i + 1L
      if (i > length(args) || startsWith(args[i], "--")) stop("Missing value for --", name)
      opts[[name]] <- args[i]
    }
    i <- i + 1L
  }
  if (is.null(opts$out) || !nzchar(opts$out)) stop("--out is required.")
  if (!nzchar(opts$data_dir)) stop("Supply --data-dir or CORTREE_DNA_DATA.")
  stopifnot(opts$dataset %in% c("REST", "NRF1"), opts$init %in% c("all", "rowsum", "pwm", "tss_dist"))
  opts$out <- absolute_output(opts$out)
  checkout <- normalizePath(file.path(dirname(script_path()), "..", ".."), winslash = "/")
  if (identical(opts$out, checkout) || startsWith(opts$out, paste0(checkout, "/")))
    stop("Choose --out outside the package checkout.")
  opts$data_dir <- normalizePath(opts$data_dir, mustWork = TRUE)
  opts
}

prior_catalog <- function(n, dataset) {
  make <- function(rate = 0, upper = Inf, a = 0, hierarchy = FALSE, independent = FALSE)
    list(ghs_diag_rate = rate, ghs_diag_upper = upper, ghs_det_df = a,
      ghs_scale_hierarchy = hierarchy, ghs_scale_shape = 3,
      ghs_scale_rate_shape = 4, ghs_scale_rate_rate = 2, all_ind = independent)
  choices <- list(flat = make(), uniform = make(upper = 100),
    exponential = make(rate = log(2)), indtree = make(independent = TRUE))
  for (r in c(0.25, 1, 4)) choices[[paste0("bounded_r", r)]] <- make(upper = 100, a = r * n / 5)
  for (s0 in c(0.1, 0.3, 1)) {
    suffix <- c("0.1" = "s0p1", "0.3" = "s0p3", "1" = "s1p0")[[as.character(s0)]]
    choices[[paste0("soft_", suffix)]] <- make(rate = 2 * s0 * s0, a = 4)
    choices[[paste0("hier_", suffix)]] <- make(rate = 2 * s0 * s0, a = 4, hierarchy = TRUE)
    if (dataset == "REST" || s0 == 0.1) {
      a <- n - n / 5
      choices[[paste0("strong_soft_", suffix)]] <- make(rate = a * s0 * s0 / 2, a = a)
    }
  }
  choices
}

# Exact Dahl loss via contingency tables avoids an N by N allocation.
dahl_partition <- function(Z) {
  draws <- ncol(Z)
  codes <- lapply(seq_len(draws), function(i) match(Z[, i], unique(Z[, i])))
  K <- vapply(codes, max, integer(1L)); overlap <- matrix(0, draws, draws)
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

run_one <- function(directory, config, execute) {
  if (file.exists(file.path(directory, "COMPLETE"))) {
    stopifnot(identical(readRDS(file.path(directory, "metadata.rds"))$config, config))
    hashes <- readRDS(file.path(directory, "artifact_md5.rds"))
    stopifnot(identical(unname(tools::md5sum(file.path(directory, names(hashes)))), unname(hashes)))
    message("Verified existing result: ", basename(directory))
    return(read.csv(file.path(directory, "metrics.csv"), stringsAsFactors = FALSE))
  }
  if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
    stop("Incomplete output exists; inspect it and use a fresh --out: ", directory)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(directory, ".running")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Output is already running: ", directory)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  old <- getwd(); setwd(directory); on.exit(setwd(old), add = TRUE)
  # A baseline package that opens a graphics device must write outside checkout.
  options(device = function(...) grDevices::pdf(file = file.path(directory, "Rplots.pdf"), ...))
  warnings_seen <- character()
  tryCatch({
    result <- withCallingHandlers(execute(directory), warning = function(w) {
      warnings_seen <<- unique(c(warnings_seen, conditionMessage(w)))
    })
    result$metrics$warnings <- paste(warnings_seen, collapse = " | ")
    write.csv(result$metrics, file.path(directory, "metrics.csv"), row.names = FALSE)
    saveRDS(list(config = config, prior = result$prior, completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
      file.path(directory, "metadata.rds"))
    writeLines(warnings_seen, file.path(directory, "warnings.txt"))
    writeLines(capture.output(sessionInfo()), file.path(directory, "sessionInfo.txt"))
    files <- list.files(directory, full.names = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    hashes <- tools::md5sum(files); names(hashes) <- basename(files)
    saveRDS(hashes, file.path(directory, "artifact_md5.rds"))
    writeLines(format(Sys.time(), tz = "UTC", usetz = TRUE), file.path(directory, "COMPLETE"))
    result$metrics
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(directory, "FAILED"))
    stop(e)
  })
}

main <- function() {
  Sys.setenv(TZ = "UTC")
  opts <- read_options()
  if (!is.null(opts$library)) .libPaths(c(normalizePath(opts$library, mustWork = TRUE), .libPaths()))
  if (!requireNamespace("CorTree", quietly = TRUE)) stop("Install CorTree first; see DNASE.md.")
  if (!is.null(opts$library) && normalizePath(find.package("CorTree")) !=
      normalizePath(file.path(opts$library, "CorTree"))) stop("CorTree was not loaded from --library.")
  # Match the R 4.4.3 archived entry state. Quantile initialization draws no RNG.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  if (getRversion() != "4.4.3") message("Reference R version is 4.4.3; exact random replay across versions is not promised.")
  required_args <- c("ghs_diag_rate", "ghs_diag_upper", "ghs_det_df", "ghs_scale_hierarchy")
  if (!all(required_args %in% names(formals(CorTree::CorTree_sampler))))
    stop("This runner requires the updated CorTree precision-prior API.")
  files <- file.path(opts$data_dir, paste0(opts$dataset,
    c(".K562.DNase.counts.mat.rds", ".K562.sites.chip.labels.tss.dist.rds")))
  if (!all(file.exists(files))) stop("Missing required input: ", paste(files[!file.exists(files)], collapse = ", "))
  raw <- readRDS(files[1L]); sites <- as.data.frame(readRDS(files[2L]))
  stopifnot(is.matrix(raw), is.numeric(raw), nrow(raw) == nrow(sites), ncol(raw) %% 2L == 0L,
    all(c("pwm.score", "strand", "tss.dist", "chip_label", "chip") %in% names(sites)),
    all(is.finite(raw)), all(raw >= 0), all(raw == floor(raw)))
  p <- ncol(raw) %/% 2L
  merged <- raw[, seq_len(p), drop = FALSE] + raw[, p + seq_len(p), drop = FALSE]
  keep <- which(rowSums(merged) >= 50 & sites$pwm.score >= 13 & sites$strand == "+")
  X <- as.matrix(merged[keep, , drop = FALSE]); sites <- sites[keep, , drop = FALSE]
  expected_n <- if (opts$dataset == "REST") 4551L else 8148L
  stopifnot(nrow(X) == expected_n, ncol(X) == if (opts$dataset == "REST") 220L else 211L,
    all(is.finite(sites$pwm.score)), all(is.finite(sites$tss.dist)),
    !anyNA(sites$chip_label), all(sites$chip_label %in% c(0, 1)))
  full_inits <- list(rowsum = CorTree::quantile_groups(rowSums(X), 5L) - 1L,
    pwm = CorTree::quantile_groups(sites$pwm.score, 5L) - 1L,
    tss_dist = CorTree::quantile_groups(sites$tss.dist, 5L) - 1L)
  take <- if (opts$smoke) seq_len(min(64L, nrow(X))) else seq_len(nrow(X))
  X <- X[take, , drop = FALSE]; sites <- sites[take, , drop = FALSE]; keep <- keep[take]
  truth <- as.integer(sites$chip_label)
  rm(raw, merged); invisible(gc())
  catalog <- prior_catalog(expected_n, opts$dataset)
  selected <- c("flat", "uniform", "soft_s0p1", "strong_soft_s0p1", "indtree")
  priors <- switch(opts$prior, selected = selected, all = names(catalog), baselines = character(), opts$prior)
  if (any(!priors %in% names(catalog))) stop("Unknown/unrun prior for ", opts$dataset,
    "; choose selected, all, baselines, or: ", paste(names(catalog), collapse = ", "))
  inits <- if (opts$init == "all") names(full_inits) else opts$init
  package_dir <- find.package("CorTree")
  package_files <- list.files(package_dir, recursive = TRUE, full.names = TRUE)
  package_md5 <- tools::md5sum(package_files)
  names(package_md5) <- substring(package_files, nchar(package_dir) + 2L)
  base_config <- list(dataset = opts$dataset, full_filtered_n = expected_n, n = nrow(X), p = p,
    input_files = files, input_md5 = tools::md5sum(files), retained_input_rows = keep,
    smoke = opts$smoke, count_threshold = 50, pwm_threshold = 13, strand = "+",
    seed = 2025L, RNGkind = RNGkind(), R_version = R.version.string,
    package_version = as.character(utils::packageVersion("CorTree")), package_md5 = package_md5,
    runner_md5 = unname(tools::md5sum(script_path())))
  suite <- file.path(opts$out, if (opts$smoke) "smoke" else "production", opts$dataset)
  dir.create(suite, recursive = TRUE, showWarnings = FALSE)
  total_iter <- if (opts$smoke) 6L else 150L; burnin <- if (opts$smoke) 2L else 100L
  summaries <- list()
  for (prior_name in priors) for (init_name in inits) {
    prior_args <- catalog[[prior_name]]
    initial <- as.integer(full_inits[[init_name]][take])
    stopifnot(length(initial) == nrow(X), !anyNA(initial), all(initial %in% 0:4))
    args <- c(list(X = X, init_Z = initial, n_clus = 5L, tree_depth = 9L,
      cutoff_layer = if (opts$dataset == "REST") 3L else 4L,
      total_iter = total_iter, burnin = burnin, warm_start = 0L,
      c_sigma2_vec = 10, sigma_mu2 = 0.1, cov_interval = 3L,
      save_phi_trace = FALSE, save_cluster_cor_trace = FALSE), prior_args)
    config <- c(base_config, list(prior = prior_name, initialization = init_name,
      initial_labels = initial, sampler_args = args[setdiff(names(args), c("X", "init_Z"))]))
    directory <- file.path(suite, paste(prior_name, init_name, sep = "_"))
    summaries[[length(summaries) + 1L]] <- run_one(directory, config, function(directory) {
      saveRDS(full_inits[[init_name]], file.path(directory, "full_initial_labels.rds"))
      saveRDS(initial, file.path(directory, "initial_labels.rds"))
      set.seed(2025L)
      saveRDS(.Random.seed, file.path(directory, "rng_before_sampler.rds"))
      fit <- do.call(CorTree::CorTree_sampler, args)
      saveRDS(.Random.seed, file.path(directory, "rng_after_sampler.rds"))
      Z <- as.matrix(fit$mcmc$Z)
      stopifnot(identical(dim(Z), c(nrow(X), total_iter - burnin)), all(Z %in% 0:4),
        all(is.finite(fit$mcmc$loglik)), all(is.finite(fit$mcmc$pi)),
        max(abs(colSums(fit$mcmc$pi) - 1)) < 1e-8)
      dahl <- dahl_partition(Z)
      canonical <- vapply(seq_len(ncol(Z)), function(i) match(Z[, i], unique(Z[, i])), integer(nrow(Z)))
      guard <- fit$mcmc$covariance_safeguards
      if (prior_name != "flat" && !prior_args$all_ind) stopifnot(isTRUE(fit$mcmc$ghs_prior$proper),
        !guard$enabled, sum(guard$regularizations) == 0, sum(guard$updates_skipped) == 0)
      if (!prior_args$all_ind) stopifnot(
        isTRUE(all.equal(fit$mcmc$ghs_prior$diag_rate, prior_args$ghs_diag_rate)),
        isTRUE(all.equal(fit$mcmc$ghs_prior$diag_upper, prior_args$ghs_diag_upper)),
        isTRUE(all.equal(fit$mcmc$ghs_prior$det_df, prior_args$ghs_det_df)))
      precision <- list()
      if (!prior_args$all_ind) for (i in seq_along(fit$mcmc$Sigma_inv)) for (k in 1:5) {
        omega <- fit$mcmc$Sigma_inv[[i]][, , k]
        ev <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
        stopifnot(all(is.finite(omega)), min(ev) > 0,
          max(abs(omega - t(omega))) < 1e-8 * max(1, max(abs(omega))))
        if (is.finite(prior_args$ghs_diag_upper)) stopifnot(all(diag(omega) < prior_args$ghs_diag_upper))
        precision[[length(precision) + 1L]] <- data.frame(draw = i, component = k - 1L,
          occupancy = sum(Z[, i] == k - 1L), min_eigenvalue = min(ev), max_diagonal = max(diag(omega)))
      }
      if (length(precision)) write.csv(do.call(rbind, precision), file.path(directory, "precision_diagnostics.csv"), row.names = FALSE)
      if (prior_args$ghs_scale_hierarchy) stopifnot(isTRUE(fit$mcmc$ghs_scale$active),
        all(is.finite(fit$mcmc$ghs_scale$t)), all(fit$mcmc$ghs_scale$t > 0),
        all(is.finite(fit$mcmc$ghs_scale$b)), all(fit$mcmc$ghs_scale$b > 0))
      assignments <- data.frame(input_row = keep, chip_label = truth, init = initial, dahl = dahl$Z,
        pwm_score = sites$pwm.score, tss_dist = sites$tss.dist, total_count = rowSums(X))
      write.csv(assignments, file.path(directory, "assignments.csv"), row.names = FALSE)
      saveRDS(list(fit = fit, dahl = dahl), file.path(directory, "fit.rds"))
      list(prior = fit$mcmc$ghs_prior, metrics = data.frame(dataset = opts$dataset, prior = prior_name,
        initialization = init_name, smoke = opts$smoke, n = nrow(X), total_iter = total_iter, burnin = burnin,
        ARI = CorTree::adjusted_rand_index(truth, dahl$Z), initial_ARI = CorTree::adjusted_rand_index(truth, initial),
        n_clusters = length(unique(dahl$Z)), retained_unique_partitions = nrow(unique(t(canonical))),
        elapsed_sec = as.numeric(fit$elapsed), regularizations = sum(guard$regularizations),
        updates_skipped = sum(guard$updates_skipped)))
    })
  }
  if (opts$prior == "baselines") {
    for (pkg in c("cluster", "CENTIPEDE", "DirichletMultinomial"))
      if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing optional baseline package: ", pkg)
    for (method in c("K-means", "PAM", "CENTIPEDE", "DMM")) {
      config <- c(base_config, list(method = method, K = 2L,
        covariates = if (method == "CENTIPEDE") "PWM with intercept" else "none"))
      summaries[[length(summaries) + 1L]] <- run_one(file.path(suite, paste0("baseline_", method)), config, function(directory) {
        set.seed(2025L)
        saveRDS(.Random.seed, file.path(directory, "rng_before_sampler.rds"))
        elapsed <- system.time({
          fit <- switch(method,
            "K-means" = stats::kmeans(scale(X), centers = 2L, nstart = 25L),
            PAM = cluster::pam(X, k = 2L),
            CENTIPEDE = CENTIPEDE::fitCentipede(Xlist = list(DNase = X), Y = cbind(1, sites$pwm.score)),
            DMM = DirichletMultinomial::dmn(X, k = 2L, seed = 2025L))
        })[["elapsed"]]
        if (method == "DMM") {
          posterior <- as.matrix(DirichletMultinomial::mixture(fit))
          if (nrow(posterior) != nrow(X) && ncol(posterior) == nrow(X)) posterior <- t(posterior)
          stopifnot(identical(dim(posterior), c(nrow(X), 2L)), all(is.finite(posterior)))
          z <- max.col(posterior, ties.method = "first") - 1L
        } else z <- switch(method, "K-means" = fit$cluster, PAM = fit$clustering,
          CENTIPEDE = as.integer(fit$PostPr > 0.5))
        stopifnot(length(z) == nrow(X), !anyNA(z))
        saveRDS(.Random.seed, file.path(directory, "rng_after_sampler.rds"))
        if (method == "PAM") fit$diss <- NULL
        saveRDS(fit, file.path(directory, "fit.rds"))
        write.csv(data.frame(input_row = keep, chip_label = truth, cluster = z),
          file.path(directory, "assignments.csv"), row.names = FALSE)
        list(prior = NULL, metrics = data.frame(dataset = opts$dataset, prior = method,
          initialization = "none", smoke = opts$smoke, n = nrow(X),
          ARI = CorTree::adjusted_rand_index(truth, z), n_clusters = length(unique(z)), elapsed_sec = elapsed))
      })
    }
  }
  write.csv(do.call(rbind, summaries),
    file.path(suite, paste0("summary_", opts$prior, "_", opts$init, ".csv")), row.names = FALSE)
  cat("Completed requested fits under ", suite, "\n", sep = "")
}

main()
