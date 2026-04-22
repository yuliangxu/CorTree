resolve_data_dir <- function(project_dir, file_names) {
  candidates <- unique(Filter(
    nzchar,
    c(
      Sys.getenv("CORTREE_DNA_DATA", unset = ""),
      file.path(project_dir, "data"),
      "/hpc/group/mastatlab/yx306/DNA_data",
      "/Users/yuliangxu/Library/Mobile Documents/com~apple~CloudDocs/MyDrive/Research/Tree/code/AutoTree/data"
    )
  ))
  for (dir_path in candidates) {
    if (all(file.exists(file.path(dir_path, file_names)))) {
      return(dir_path)
    }
  }
  stop(
    "Could not resolve a data directory containing: ",
    paste(file_names, collapse = ", ")
  )
}

quantile_groups <- function(x, k) {
  qs <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  if (anyDuplicated(qs) > 0L) {
    qs <- unique(qs)
    if (length(qs) < 2L) {
      return(rep.int(1L, length(x)))
    }
    return(as.integer(cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE)))
  }
  as.integer(cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE))
}

postburnin_cols <- function(n_saved, burnin) {
  start <- min(as.integer(burnin) + 1L, n_saved)
  seq.int(start, n_saved)
}

hard_cluster_from_chain <- function(Z_chain, n_clus, burnin) {
  keep_cols <- postburnin_cols(ncol(Z_chain), burnin)
  Z_keep <- Z_chain[, keep_cols, drop = FALSE]
  Z_hat <- as.integer(round(rowMeans(Z_keep)))
  pmax(0L, pmin(as.integer(n_clus) - 1L, Z_hat))
}

extract_pi_trace <- function(fit, burnin) {
  pi_chain <- fit$mcmc$pi
  keep_cols <- postburnin_cols(ncol(pi_chain), burnin)
  pi_chain[, keep_cols, drop = FALSE]
}

adjusted_rand_index <- function(labels_true, labels_pred) {
  if (length(labels_true) != length(labels_pred)) {
    stop("labels_true and labels_pred must have the same length")
  }

  ok <- !(is.na(labels_true) | is.na(labels_pred))
  labels_true <- labels_true[ok]
  labels_pred <- labels_pred[ok]
  n <- length(labels_true)

  if (n < 2L) {
    return(1)
  }

  contingency <- table(labels_true, labels_pred)
  choose2 <- function(x) x * (x - 1) / 2

  sum_ij <- sum(choose2(contingency))
  sum_i <- sum(choose2(rowSums(contingency)))
  sum_j <- sum(choose2(colSums(contingency)))
  total_pairs <- choose2(n)

  if (total_pairs == 0) {
    return(1)
  }

  expected_index <- (sum_i * sum_j) / total_pairs
  max_index <- 0.5 * (sum_i + sum_j)
  denom <- max_index - expected_index

  if (isTRUE(all.equal(denom, 0))) {
    return(1)
  }

  (sum_ij - expected_index) / denom
}

extract_task_id <- function(path) {
  base_name <- basename(path)
  task_match <- regmatches(base_name, regexec("(^|_)task([0-9]+)_", base_name))[[1]]
  if (length(task_match) < 3L) {
    return(NA_integer_)
  }
  as.integer(task_match[3])
}

extract_elapsed <- function(bundle, method_key, init_mode = NULL) {
  timing <- bundle$timing
  if (is.null(timing)) {
    if (identical(method_key, "CorTree") && !is.null(bundle$fit$cortree$elapsed)) {
      return(as.numeric(bundle$fit$cortree$elapsed))
    }
    if (identical(method_key, "IndTree") && !is.null(bundle$fit$indtree$elapsed)) {
      return(as.numeric(bundle$fit$indtree$elapsed))
    }
    return(NA_real_)
  }

  if (identical(method_key, "Z_init")) {
    if (is.null(init_mode) || is.null(timing$init[[init_mode]])) {
      return(NA_real_)
    }
    return(as.numeric(timing$init[[init_mode]]))
  }

  if (identical(method_key, "CorTree")) {
    return(as.numeric(if (!is.null(timing$cortree)) timing$cortree else bundle$fit$cortree$elapsed))
  }
  if (identical(method_key, "IndTree")) {
    return(as.numeric(if (!is.null(timing$indtree)) timing$indtree else bundle$fit$indtree$elapsed))
  }
  if (identical(method_key, "K-means") && !is.null(timing$kmeans)) {
    return(as.numeric(timing$kmeans))
  }
  if (identical(method_key, "PAM") && !is.null(timing$pam)) {
    return(as.numeric(timing$pam))
  }
  if (identical(method_key, "CENTIPEDE") && !is.null(timing$centipede)) {
    return(as.numeric(timing$centipede))
  }

  NA_real_
}

extract_sensitivity_mode <- function(path) {
  base_name <- basename(path)
  mode_match <- regmatches(
    base_name,
    regexec("_task[0-9]+_(.*)_result_bundle\\.rds$", base_name)
  )[[1]]
  if (length(mode_match) < 2L) {
    return(NA_character_)
  }
  mode_match[2]
}

read_sensitivity_bundle_records <- function(dataset, out_dir, sensitivity_label_map) {
  bundle_paths <- list.files(
    out_dir,
    pattern = sprintf("^%s_sensi_task[0-9]+_.*_result_bundle\\.rds$", dataset),
    full.names = TRUE
  )

  if (length(bundle_paths) == 0L) {
    return(NULL)
  }

  bundle_paths <- bundle_paths[order(vapply(bundle_paths, extract_task_id, integer(1)))]

  records <- lapply(bundle_paths, function(path) {
    bundle <- readRDS(path)
    tab <- as.data.frame(bundle$ari_table, stringsAsFactors = FALSE)
    required_cols <- c("Method", "ARI")
    missing_cols <- setdiff(required_cols, names(tab))
    if (length(missing_cols) > 0L) {
      stop(
        dataset,
        " bundle ari_table is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        " in ",
        basename(path)
      )
    }

    sensitivity_mode <- if (!is.null(bundle$sensitivity_mode)) {
      as.character(bundle$sensitivity_mode)
    } else {
      extract_sensitivity_mode(path)
    }
    sensitivity_label <- if (!is.null(bundle$sensitivity_label)) {
      as.character(bundle$sensitivity_label)
    } else {
      unname(sensitivity_label_map[[sensitivity_mode]])
    }

    data.frame(
      Dataset = dataset,
      TaskID = rep.int(extract_task_id(path), nrow(tab)),
      SensitivityMode = rep.int(sensitivity_mode, nrow(tab)),
      SensitivityLabel = rep.int(sensitivity_label, nrow(tab)),
      Method = as.character(tab$Method),
      MethodKey = ifelse(grepl("^Z_init", tab$Method), "Z_init", as.character(tab$Method)),
      ARI = as.numeric(tab$ARI),
      SourceType = "result_bundle",
      SourceFile = rep.int(basename(path), nrow(tab)),
      BundlePath = rep.int(path, nrow(tab)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind.data.frame, c(records, list(stringsAsFactors = FALSE)))
}

read_sensitivity_rout_records <- function(dataset, out_dir, sensitivity_label_map) {
  rout_paths <- list.files(
    out_dir,
    pattern = sprintf("^%s_sensi[0-9]+\\.Rout$", dataset),
    full.names = TRUE
  )

  if (length(rout_paths) == 0L) {
    return(NULL)
  }

  rout_paths <- rout_paths[order(as.integer(gsub("^.*sensi([0-9]+)\\.Rout$", "\\1", rout_paths)))]

  records <- lapply(rout_paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    sensitivity_line <- grep("^Sensitivity:", lines, value = TRUE)
    if (length(sensitivity_line) == 0L) {
      stop("Could not find sensitivity label in ", basename(path))
    }
    sensitivity_label <- trimws(sub("^Sensitivity:\\s*", "", tail(sensitivity_line, 1L)))
    sensitivity_mode_match <- names(sensitivity_label_map)[match(
      sensitivity_label,
      unname(sensitivity_label_map)
    )]
    sensitivity_mode <- if (length(sensitivity_mode_match) == 0L || is.na(sensitivity_mode_match)) {
      NA_character_
    } else {
      sensitivity_mode_match
    }

    method_idx <- grep("^\\s*Method\\s+ARI\\s*$", lines)
    if (length(method_idx) == 0L) {
      stop("Could not find ARI table header in ", basename(path))
    }
    table_lines <- lines[seq.int(tail(method_idx, 1L) + 1L, min(length(lines), tail(method_idx, 1L) + 2L))]

    parsed_rows <- lapply(table_lines, function(line) {
      cleaned <- trimws(line)
      parts <- strsplit(cleaned, "\\s+", perl = TRUE)[[1]]
      ari_value <- suppressWarnings(as.numeric(tail(parts, 1L)))
      if (length(parts) < 3L || is.na(ari_value)) {
        stop("Could not parse ARI row in ", basename(path), ": ", line)
      }
      data.frame(
        Method = paste(parts[2:(length(parts) - 1L)], collapse = " "),
        ARI = ari_value,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    tab <- do.call(rbind.data.frame, c(parsed_rows, list(stringsAsFactors = FALSE)))

    data.frame(
      Dataset = dataset,
      TaskID = rep.int(as.integer(gsub("^.*sensi([0-9]+)\\.Rout$", "\\1", basename(path))), nrow(tab)),
      SensitivityMode = rep.int(sensitivity_mode, nrow(tab)),
      SensitivityLabel = rep.int(sensitivity_label, nrow(tab)),
      Method = as.character(tab$Method),
      MethodKey = ifelse(grepl("^Z_init", tab$Method), "Z_init", as.character(tab$Method)),
      ARI = as.numeric(tab$ARI),
      SourceType = "rout",
      SourceFile = rep.int(basename(path), nrow(tab)),
      BundlePath = rep.int(NA_character_, nrow(tab)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  do.call(rbind.data.frame, c(records, list(stringsAsFactors = FALSE)))
}

collect_sensitivity_dataset_results <- function(dataset, out_dir, fallback_dir, sensitivity_label_map) {
  bundle_records <- read_sensitivity_bundle_records(dataset, out_dir, sensitivity_label_map)
  if (!is.null(bundle_records) && nrow(bundle_records) > 0L) {
    return(bundle_records[order(bundle_records$TaskID, bundle_records$MethodKey), , drop = FALSE])
  }

  rout_records <- read_sensitivity_rout_records(dataset, fallback_dir, sensitivity_label_map)
  if (!is.null(rout_records) && nrow(rout_records) > 0L) {
    return(rout_records[order(rout_records$TaskID, rout_records$MethodKey), , drop = FALSE])
  }

  stop("No sensitivity results found for ", dataset, " in ", out_dir, " or ", fallback_dir)
}

filter_sensitivity_records <- function(records) {
  if (nrow(records) == 0L) {
    return(records)
  }

  bundle_rows <- !is.na(records$BundlePath) & nzchar(records$BundlePath)
  if (!any(bundle_rows)) {
    return(records)
  }

  bundle_keys <- unique(records[bundle_rows, c("Dataset", "TaskID", "BundlePath"), drop = FALSE])
  if (nrow(bundle_keys) == 0L) {
    return(records)
  }

  keep_paths <- unlist(lapply(split(bundle_keys, paste(bundle_keys$Dataset, bundle_keys$TaskID, sep = "::")), function(key_df) {
    if (nrow(key_df) <= 1L) {
      return(key_df$BundlePath[1])
    }

    path_scores <- vapply(key_df$BundlePath, function(path) {
      cor_row <- records[
        records$BundlePath == path & records$MethodKey == "CorTree",
        ,
        drop = FALSE
      ]
      if (nrow(cor_row) == 0L) {
        return(NA_real_)
      }
      cor_row$ARI[1]
    }, numeric(1))

    key_df$BundlePath[which.max(path_scores)]
  }), use.names = FALSE)

  keep_rows <- !bundle_rows | records$BundlePath %in% keep_paths
  records[keep_rows, , drop = FALSE]
}

build_sensitivity_task_summary <- function(records) {
  task_info <- unique(records[, c("Dataset", "TaskID", "SensitivityMode", "SensitivityLabel", "BundlePath"), drop = FALSE])
  task_info <- task_info[order(task_info$TaskID, task_info$SensitivityMode), , drop = FALSE]
  row.names(task_info) <- NULL

  key_cols <- c("Dataset", "TaskID", "SensitivityMode", "SensitivityLabel", "BundlePath")
  z_init <- records[records$MethodKey == "Z_init", c(key_cols, "ARI"), drop = FALSE]
  names(z_init)[length(names(z_init))] <- "Z_init_ARI"
  cortree <- records[records$MethodKey == "CorTree", c(key_cols, "ARI"), drop = FALSE]
  names(cortree)[length(names(cortree))] <- "CorTree_ARI"

  out <- merge(task_info, z_init, by = key_cols, all.x = TRUE, sort = FALSE)
  out <- merge(out, cortree, by = key_cols, all.x = TRUE, sort = FALSE)
  out$ARI_Improvement <- out$CorTree_ARI - out$Z_init_ARI
  out <- out[order(out$TaskID, out$SensitivityMode), , drop = FALSE]
  row.names(out) <- NULL
  out
}

build_sensitivity_method_wide <- function(records) {
  sensitivity_order <- unique(records$SensitivityMode[order(records$TaskID)])
  ari_matrix <- with(
    records,
    tapply(ARI, list(MethodKey, SensitivityMode), identity)
  )
  out <- data.frame(
    MethodKey = row.names(ari_matrix),
    as.data.frame.matrix(ari_matrix),
    row.names = NULL,
    check.names = FALSE
  )
  wide_cols <- intersect(sensitivity_order, names(out))
  out <- out[, c("MethodKey", wide_cols), drop = FALSE]
  out$MeanARI <- rowMeans(as.matrix(out[, wide_cols, drop = FALSE]), na.rm = TRUE)
  out <- out[order(-out$MeanARI), , drop = FALSE]
  row.names(out) <- NULL
  out
}

format_sensitivity_dataset_summary <- function(task_summary) {
  best_idx <- which.max(task_summary$CorTree_ARI)
  worst_idx <- which.min(task_summary$CorTree_ARI)
  mean_improvement <- mean(task_summary$ARI_Improvement, na.rm = TRUE)

  c(
    sprintf("Dataset: %s", unique(task_summary$Dataset)),
    sprintf(
      "Baseline Z_init ARI is constant across sensitivity settings: %.4f",
      unique(task_summary$Z_init_ARI)
    ),
    sprintf(
      "CorTree ARI ranges from %.4f (%s) to %.4f (%s).",
      task_summary$CorTree_ARI[worst_idx],
      task_summary$SensitivityLabel[worst_idx],
      task_summary$CorTree_ARI[best_idx],
      task_summary$SensitivityLabel[best_idx]
    ),
    sprintf(
      "Mean CorTree improvement over Z_init across settings: %.4f ARI.",
      mean_improvement
    ),
    ""
  )
}

reconstruct_bundle_inputs <- function(bundle) {
  data_dir <- bundle$input$data_dir
  file_name <- bundle$input$file_name
  label_name <- bundle$input$label_name

  dnase_count_matrix <- readRDS(file.path(data_dir, file_name))
  sites_chip_labels <- readRDS(file.path(data_dir, label_name))
  dnase_count_matrix <- dnase_count_matrix[, 1:(ncol(dnase_count_matrix) / 2)] +
    dnase_count_matrix[, (ncol(dnase_count_matrix) / 2 + 1):ncol(dnase_count_matrix)]

  idx_include <- bundle$input$idx_include
  if (is.null(idx_include)) {
    threshold <- if (!is.null(bundle$input$threshold)) bundle$input$threshold else 50
    pwm_score_cutoff <- if (!is.null(bundle$input$pwm_score_cutoff)) bundle$input$pwm_score_cutoff else 13
    idx_include <- which(
      rowSums(dnase_count_matrix) >= threshold &
        sites_chip_labels$pwm.score >= pwm_score_cutoff &
        sites_chip_labels$strand == "+"
    )
  }

  list(
    X = as.matrix(dnase_count_matrix[idx_include, , drop = FALSE]),
    chip_labels = as.numeric(sites_chip_labels$chip_label[idx_include]),
    chip_data = sites_chip_labels$chip[idx_include]
  )
}

load_bundle_plot_data <- function(path) {
  bundle <- readRDS(path)
  reconstructed <- reconstruct_bundle_inputs(bundle)
  out <- list(
    X = if (!is.null(bundle$data$X)) bundle$data$X else reconstructed$X,
    chip_labels = if (!is.null(bundle$data$chip_labels)) bundle$data$chip_labels else reconstructed$chip_labels,
    chip_data = if (!is.null(bundle$data$chip_data)) bundle$data$chip_data else reconstructed$chip_data,
    cortree = bundle$cluster_assignments$cortree,
    init_mode = if (!is.null(bundle$init_mode)) bundle$init_mode else NA_character_,
    sensitivity_label = if (!is.null(bundle$sensitivity_label)) bundle$sensitivity_label else NA_character_
  )
  rm(bundle, reconstructed)
  gc()
  out
}

resize_plot_text <- function(plot_obj) {
  plot_obj +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 17, face = "plain", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.size = grid::unit(0.6, "cm"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

rename_plot_title <- function(plot_obj, title_text) {
  plot_obj + ggplot2::labs(title = title_text)
}

select_final_bundle_path <- function(dataset, out_dir) {
  preferred_path <- file.path(out_dir, sprintf("%s_task02_pwm_quantile_result_bundle.rds", dataset))
  if (file.exists(preferred_path)) {
    return(preferred_path)
  }

  candidate_paths <- list.files(
    out_dir,
    pattern = sprintf("^%s_task[0-9]+_.*_result_bundle\\.rds$", dataset),
    full.names = TRUE
  )
  if (length(candidate_paths) == 0L) {
    stop("No final analysis result bundles found for ", dataset, " in ", out_dir)
  }

  cor_ari <- vapply(candidate_paths, function(path) {
    bundle <- readRDS(path)
    tab <- as.data.frame(bundle$ari_table, stringsAsFactors = FALSE)
    idx <- which(tab$Method == "CorTree")
    if (length(idx) == 0L) {
      return(NA_real_)
    }
    as.numeric(tab$ARI[idx[1]])
  }, numeric(1))

  candidate_paths[which.max(cor_ari)]
}

build_combined_sensitivity_dataset_plot <- function(dataset, records, out_dir) {
  task_summary <- build_sensitivity_task_summary(records)
  task_summary <- task_summary[order(task_summary$TaskID), , drop = FALSE]

  if (any(is.na(task_summary$BundlePath) | !nzchar(task_summary$BundlePath))) {
    stop("Combined sensitivity plots require result_bundle records for ", dataset)
  }

  selected_task_paths <- task_summary$BundlePath[match(1:3, task_summary$TaskID)]
  if (any(is.na(selected_task_paths))) {
    stop("Missing one or more sensitivity task bundles for ", dataset)
  }

  selected_task_info <- lapply(selected_task_paths, load_bundle_plot_data)
  final_bundle_path <- select_final_bundle_path(dataset, out_dir)
  final_bundle_info <- load_bundle_plot_data(final_bundle_path)
  base_info <- selected_task_info[[1L]]

  max_hist_y <- max(
    vapply(
      c(selected_task_info, list(final_bundle_info)),
      function(info) {
        h <- hist(log(info$chip_data + 1), plot = FALSE, breaks = 30)
        max(h$counts, na.rm = TRUE)
      },
      numeric(1)
    ),
    na.rm = TRUE
  )

  max_signal_y <- NULL

  cluster_titles <- c(
    "ChIP_label",
    vapply(seq_along(selected_task_info), function(i) {
      label_text <- task_summary$SensitivityLabel[task_summary$TaskID == i][1]
      if (identical(dataset, "NRF1") && identical(i, 1L)) {
        label_text <- "cutoff_layer=3"
      }
      sprintf("Task %d\n%s", i, label_text)
    }, character(1)),
    sprintf("Final result\n%s", final_bundle_info$init_mode)
  )

  cluster_mean_panels <- c(
    list(
      rename_plot_title(
        resize_plot_text(cluster_mean_gg(base_info$X, base_info$chip_labels, title = "ChIP_label", y_max = max_signal_y)),
        cluster_titles[[1]]
      )
    ),
    lapply(seq_along(selected_task_info), function(i) {
      info <- selected_task_info[[i]]
      rename_plot_title(
        resize_plot_text(cluster_mean_gg(
          info$X,
          info$cortree,
          title = paste0("Task ", i),
          y_max = max_signal_y
        )),
        cluster_titles[[i + 1L]]
      )
    }),
    list(
      rename_plot_title(
        resize_plot_text(cluster_mean_gg(
          final_bundle_info$X,
          final_bundle_info$cortree,
          title = "Final result",
          y_max = max_signal_y
        )),
        cluster_titles[[5]]
      )
    )
  )

  histogram_panels <- c(
    list(
      rename_plot_title(
        resize_plot_text(hist_by_group_gg(base_info$chip_data, base_info$chip_labels, title = "ChIP_label", y_max = max_hist_y)),
        cluster_titles[[1]]
      )
    ),
    lapply(seq_along(selected_task_info), function(i) {
      info <- selected_task_info[[i]]
      rename_plot_title(
        resize_plot_text(hist_by_group_gg(
          info$chip_data,
          info$cortree,
          title = paste0("Task ", i),
          y_max = max_hist_y
        )),
        cluster_titles[[i + 1L]]
      )
    }),
    list(
      rename_plot_title(
        resize_plot_text(hist_by_group_gg(
          final_bundle_info$chip_data,
          final_bundle_info$cortree,
          title = "Final result",
          y_max = max_hist_y
        )),
        cluster_titles[[5]]
      )
    )
  )

  patchwork::wrap_plots(c(cluster_mean_panels, histogram_panels), ncol = 5) +
    patchwork::plot_annotation(title = sprintf("%s DNase Sensitivity Summary", dataset))
}

extract_dmm_posterior <- function(dmm_fit) {
  posterior <- tryCatch(
    as.matrix(DirichletMultinomial::mixture(dmm_fit)),
    error = function(e) NULL
  )
  if (is.null(posterior)) {
    stop("Could not extract per-sample posterior probabilities from the DMM fit.")
  }
  if (nrow(posterior) < ncol(posterior)) {
    posterior <- t(posterior)
  }
  storage.mode(posterior) <- "numeric"
  posterior
}

make_sim_summary_table <- function(ari_mat, time_mat, method_names, experiment_label) {
  data.frame(
    Experiment = rep.int(experiment_label, length(method_names)),
    Method = method_names,
    Mean_ARI = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, mean, na.rm = TRUE),
    SD_ARI = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, stats::sd, na.rm = TRUE),
    Median_ARI = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, stats::median, na.rm = TRUE),
    N_ARI_eq_0 = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, function(x) sum(x == 0, na.rm = TRUE)),
    N_ARI_eq_1 = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, function(x) sum(x == 1, na.rm = TRUE)),
    Avg_time_sec = apply(time_mat[, seq_along(method_names), drop = FALSE], 2, mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

format_latex_value <- function(x, digits = 2) {
  formatC(x, digits = digits, format = "f")
}

build_sim_summary_latex_table <- function(summary_df) {
  lines <- c(
    "\\begin{tabular}{llrrrrrr}",
    "\\hline",
    "Experiment & Method & Mean ARI & S.D.\\ ARI & Median ARI & \\# of ARI=0 & \\# of ARI=1 & Avg.\\ time (sec) \\\\",
    "\\hline"
  )

  experiment_levels <- unique(summary_df$Experiment)
  for (exp_label in experiment_levels) {
    exp_df <- summary_df[summary_df$Experiment == exp_label, , drop = FALSE]
    best_mean <- max(exp_df$Mean_ARI, na.rm = TRUE)
    best_n_one <- max(exp_df$N_ARI_eq_1, na.rm = TRUE)

    for (i in seq_len(nrow(exp_df))) {
      row <- exp_df[i, , drop = FALSE]
      exp_cell <- if (i == 1L) {
        paste0("\\multirow{", nrow(exp_df), "}{*}{$m=", row$Experiment, "$}")
      } else {
        ""
      }
      mean_str <- format_latex_value(row$Mean_ARI)
      n_one_str <- as.character(row$N_ARI_eq_1)
      if (isTRUE(all.equal(row$Mean_ARI, best_mean))) {
        mean_str <- paste0("\\textbf{", mean_str, "}")
      }
      if (identical(row$N_ARI_eq_1, best_n_one)) {
        n_one_str <- paste0("\\textbf{", n_one_str, "}")
      }
      lines <- c(
        lines,
        paste0(
          exp_cell, " & ", row$Method,
          " & ", mean_str,
          " & ", format_latex_value(row$SD_ARI),
          " & ", format_latex_value(row$Median_ARI),
          " & ", row$N_ARI_eq_0,
          " & ", n_one_str,
          " & ", format_latex_value(row$Avg_time_sec),
          " \\\\"
        )
      )
    }
    lines <- c(lines, "\\hline")
  }

  c(lines, "\\end{tabular}")
}
