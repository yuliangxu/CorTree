ntile_base <- function(x, n) {
  as.integer(cut(rank(x, ties.method = "first"), breaks = n, labels = FALSE))
}

clean_covariate <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.character(x)) {
    x_trim <- trimws(x)
    missing_tokens <- c(
      "", "unknown", "no_data", "not provided", "not applicable",
      "unspecified", "missing", "missing: not provided",
      "missing: not applicable", "missing: restricted access", "na", "n/a"
    )
    x_trim[tolower(x_trim) %in% missing_tokens] <- NA_character_
    return(x_trim)
  }
  x
}

format_numeric_label <- function(x) {
  x_chr <- format(x, trim = TRUE, scientific = FALSE, digits = 15)
  x_chr <- sub("(\\.[0-9]*?)0+$", "\\1", x_chr)
  x_chr <- sub("\\.$", "", x_chr)
  x_chr
}

canonicalize_numeric_strings <- function(x) {
  x_clean <- clean_covariate(x)
  report <- data.frame()

  if (!is.character(x_clean)) {
    return(list(value = x_clean, report = report))
  }

  x_num <- suppressWarnings(as.numeric(x_clean))
  keep <- !is.na(x_clean) & !is.na(x_num)
  x_canonical <- x_clean
  x_canonical[keep] <- format_numeric_label(x_num[keep])

  if (any(keep)) {
    raw_map <- split(x_clean[keep], x_canonical[keep])
    raw_map <- Filter(function(values) length(unique(values)) > 1L, raw_map)
    if (length(raw_map) > 0L) {
      report <- do.call(
        rbind,
        lapply(names(raw_map), function(value) {
          raw_values <- sort(unique(raw_map[[value]]))
          data.frame(
            parsed_value = suppressWarnings(as.numeric(value)),
            canonical_value = value,
            raw_values = paste(raw_values, collapse = " | "),
            n_raw_values = length(raw_values),
            stringsAsFactors = FALSE
          )
        })
      )
    }
  }

  list(value = x_canonical, report = report)
}

as_numeric_if_possible <- function(x, threshold = 0.9) {
  if (is.numeric(x) || is.integer(x)) {
    return(as.numeric(x))
  }
  if (is.logical(x)) {
    return(as.numeric(x))
  }
  if (!is.character(x)) {
    return(NULL)
  }

  x_num <- suppressWarnings(as.numeric(x))
  keep <- !is.na(x)
  if (!any(keep)) {
    return(NULL)
  }
  if (mean(!is.na(x_num[keep])) >= threshold) {
    return(x_num)
  }
  NULL
}

normalize_numeric_like_covariate <- function(x, threshold = 0.9) {
  x_clean_info <- canonicalize_numeric_strings(x)
  x_clean <- x_clean_info$value
  x_num <- as_numeric_if_possible(x_clean, threshold = threshold)

  if (is.null(x_num)) {
    return(list(
      value = x_clean,
      report = x_clean_info$report
    ))
  }

  list(
    value = x_num,
    report = x_clean_info$report
  )
}

clean_covariate_dataframe <- function(covariates, numeric_threshold = 0.9) {
  cov_df <- as.data.frame(covariates, stringsAsFactors = FALSE)
  cleaned_df <- cov_df
  report_list <- list()
  report_idx <- 0L

  for (nm in names(cov_df)) {
    if (nm == "SampleID" || grepl("(^|_)ID$", nm, ignore.case = TRUE)) {
      cleaned_df[[nm]] <- clean_covariate(cov_df[[nm]])
      next
    }

    cleaned_col <- normalize_numeric_like_covariate(cov_df[[nm]], threshold = numeric_threshold)
    cleaned_df[[nm]] <- cleaned_col$value

    if (nrow(cleaned_col$report) > 0L) {
      report_idx <- report_idx + 1L
      report_list[[report_idx]] <- cbind(
        covariate = nm,
        cleaned_col$report,
        stringsAsFactors = FALSE
      )
    }
  }

  normalization_report <- if (report_idx > 0L) {
    do.call(rbind, report_list[seq_len(report_idx)])
  } else {
    data.frame(
      covariate = character(),
      parsed_value = numeric(),
      canonical_value = character(),
      raw_values = character(),
      n_raw_values = integer(),
      stringsAsFactors = FALSE
    )
  }

  rownames(cleaned_df) <- rownames(cov_df)
  rownames(normalization_report) <- NULL

  list(
    covariates = cleaned_df,
    normalization_report = normalization_report
  )
}

analyze_covariate_association <- function(cluster, covariates, sample_ids, min_obs = 30L) {
  cov_df <- as.data.frame(covariates, stringsAsFactors = FALSE)
  row_idx <- match(sample_ids, cov_df$SampleID)
  cov_df <- cov_df[row_idx, , drop = FALSE]

  cluster_factor <- factor(cluster)
  results <- vector("list", length = ncol(cov_df))
  result_idx <- 0L

  for (nm in names(cov_df)) {
    if (nm == "SampleID") {
      next
    }

    x_raw <- clean_covariate(cov_df[[nm]])
    x_num <- as_numeric_if_possible(x_raw)

    if (!is.null(x_num)) {
      keep <- !is.na(x_num) & !is.na(cluster_factor)
      n_obs <- sum(keep)
      if (n_obs < min_obs || length(unique(cluster_factor[keep])) < 2L) {
        next
      }

      kw <- suppressWarnings(kruskal.test(x_num[keep] ~ cluster_factor[keep]))
      k <- length(unique(cluster_factor[keep]))
      epsilon_sq <- if (n_obs > k) {
        max(0, unname((kw$statistic - k + 1) / (n_obs - k)))
      } else {
        NA_real_
      }

      result_idx <- result_idx + 1L
      results[[result_idx]] <- data.frame(
        covariate = nm,
        var_type = "numeric",
        n_obs = n_obs,
        n_levels = NA_integer_,
        method = "Kruskal-Wallis",
        statistic = unname(kw$statistic),
        p_value = kw$p.value,
        effect_size = epsilon_sq,
        stringsAsFactors = FALSE
      )
      next
    }

    x_chr <- as.character(x_raw)
    keep <- !is.na(x_chr) & !is.na(cluster_factor)
    n_obs <- sum(keep)
    if (n_obs < min_obs) {
      next
    }

    x_fac <- factor(x_chr[keep])
    n_levels <- nlevels(x_fac)
    if (n_levels < 2L) {
      next
    }
    if (n_levels > min(50L, floor(0.2 * n_obs))) {
      next
    }

    tab <- table(x_fac, cluster_factor[keep])
    if (min(dim(tab)) < 2L) {
      next
    }

    chi <- suppressWarnings(chisq.test(tab))
    cramer_v <- sqrt(unname(chi$statistic) / (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L)))

    result_idx <- result_idx + 1L
    results[[result_idx]] <- data.frame(
      covariate = nm,
      var_type = "categorical",
      n_obs = n_obs,
      n_levels = n_levels,
      method = "Chi-squared",
      statistic = unname(chi$statistic),
      p_value = chi$p.value,
      effect_size = cramer_v,
      stringsAsFactors = FALSE
    )
  }

  assoc_df <- do.call(rbind, results[seq_len(result_idx)])
  if (is.null(assoc_df) || nrow(assoc_df) == 0L) {
    return(data.frame())
  }

  assoc_df$p_adj <- p.adjust(assoc_df$p_value, method = "BH")
  assoc_df$rank <- rank(-assoc_df$effect_size, ties.method = "min")
  assoc_df <- assoc_df[order(-assoc_df$effect_size, assoc_df$p_value, assoc_df$covariate), ]
  rownames(assoc_df) <- NULL
  assoc_df
}

get_pi_trace_plot_data <- function(fit, sampler_args = NULL) {
  pi_chain <- NULL
  iteration_index <- integer()
  trace_source <- NULL

  if (!is.null(fit$mcmc$pi_full)) {
    pi_chain <- as.matrix(fit$mcmc$pi_full)
    iteration_index <- seq.int(0L, ncol(pi_chain) - 1L)
    trace_source <- "full"
  } else if (!is.null(fit$mcmc$pi)) {
    pi_chain <- as.matrix(fit$mcmc$pi)
    burnin <- if (!is.null(sampler_args$burnin)) as.integer(sampler_args$burnin) else NA_integer_
    total_iter <- if (!is.null(sampler_args$total_iter)) as.integer(sampler_args$total_iter) else NA_integer_
    if (!is.na(burnin) && !is.na(total_iter) && ncol(pi_chain) == (total_iter - burnin)) {
      iteration_index <- seq.int(burnin, total_iter - 1L)
      trace_source <- "post_burnin_only"
    } else {
      iteration_index <- seq.int(0L, ncol(pi_chain) - 1L)
      trace_source <- "saved_draws"
    }
  } else {
    stop("fit$mcmc does not contain a pi trace.")
  }

  storage.mode(pi_chain) <- "numeric"
  list(
    pi_chain = pi_chain,
    iteration_index = iteration_index,
    trace_source = trace_source
  )
}

build_phylo_node_plot <- function(tree, count_data, sample_id = NULL, n_label = 10L) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop(
      "Package 'ape' is required for phylogenetic tree plotting. ",
      "Install it first with install.packages('ape')."
    )
  }
  if (!inherits(tree, "phylo")) {
    stop("`tree` must inherit from class 'phylo'.")
  }
  if (is.null(colnames(count_data))) {
    stop("`count_data` must have OTU names as column names.")
  }
  if (is.null(rownames(count_data))) {
    stop("`count_data` must have sample IDs as row names.")
  }

  sample_totals <- rowSums(count_data, na.rm = TRUE)
  if (is.null(sample_id)) {
    median_rank <- order(abs(sample_totals - stats::median(sample_totals)), sample_totals)[1L]
    sample_id <- rownames(count_data)[median_rank]
  }
  if (!sample_id %in% rownames(count_data)) {
    stop("Selected sample_id was not found in rownames(count_data): ", sample_id)
  }

  tip_idx <- match(tree$tip.label, colnames(count_data))
  if (anyNA(tip_idx)) {
    missing_tip <- tree$tip.label[is.na(tip_idx)][1L]
    stop("Tree tip label was not found in count_data columns: ", missing_tip)
  }

  tree_post <- ape::reorder.phylo(tree, order = "postorder")
  n_tip <- length(tree_post$tip.label)
  n_node <- tree_post$Nnode
  total_nodes <- n_tip + n_node
  node_values <- numeric(total_nodes)
  node_values[seq_len(n_tip)] <- as.numeric(count_data[sample_id, tip_idx])

  edge_mat <- tree_post$edge
  for (edge_i in seq.int(nrow(edge_mat), 1L)) {
    parent_i <- edge_mat[edge_i, 1L]
    child_i <- edge_mat[edge_i, 2L]
    node_values[parent_i] <- node_values[parent_i] + node_values[child_i]
  }

  root_node <- setdiff(edge_mat[, 1L], edge_mat[, 2L])[1L]
  node_values[root_node] <- sum(node_values[seq_len(n_tip)])
  internal_node_ids <- seq.int(n_tip + 1L, total_nodes)
  plot_order_values <- c(node_values[seq_len(n_tip)], node_values[internal_node_ids])

  nonzero_values <- plot_order_values[plot_order_values > 0]
  palette_fn <- grDevices::colorRampPalette(c("#6FA8DC", "#B7D97A", "#F6E945"))
  if (length(nonzero_values) > 1L) {
    breaks <- unique(stats::quantile(nonzero_values, probs = seq(0, 1, length.out = 101), na.rm = TRUE))
  } else {
    breaks <- c(0, max(1, nonzero_values))
  }
  color_bins <- cut(
    plot_order_values,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  color_bins[is.na(color_bins)] <- 1L
  node_palette <- palette_fn(max(color_bins))
  node_cols <- node_palette[color_bins]
  node_cols[plot_order_values <= 0] <- "#D9D9D9"

  node_cex <- scales::rescale(log1p(plot_order_values), to = c(0.4, 1.2))
  node_cex[!is.finite(node_cex)] <- 0.4

  internal_values <- node_values[(n_tip + 1L):total_nodes]
  label_n <- min(as.integer(n_label), length(internal_values))
  internal_rank <- order(internal_values, decreasing = TRUE)[seq_len(label_n)]
  internal_rank <- internal_rank[internal_values[internal_rank] > 0]
  internal_ids <- n_tip + internal_rank
  internal_labels <- seq_along(internal_ids)

  list(
    tree = tree_post,
    sample_id = sample_id,
    sample_total = unname(sample_totals[sample_id]),
    node_values = node_values,
    plot_values = plot_order_values,
    tip_cols = node_cols[seq_len(n_tip)],
    internal_cols = node_cols[internal_node_ids],
    tip_cex = node_cex[seq_len(n_tip)],
    internal_cex = node_cex[internal_node_ids],
    internal_node_ids = internal_node_ids,
    internal_ids = internal_ids,
    internal_labels = internal_labels
  )
}

save_phylo_node_plot <- function(plot_obj, pdf_path, png_path, structure_only = FALSE) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop(
      "Package 'ape' is required for phylogenetic tree plotting. ",
      "Install it first with install.packages('ape')."
    )
  }
  draw_once <- function() {
    ape::plot.phylo(
      plot_obj$tree,
      type = "phylogram",
      show.tip.label = FALSE,
      no.margin = TRUE,
      cex = 0.35,
      edge.color = "grey40",
      edge.width = 0.7
    )
    tip_bg <- if (isTRUE(structure_only)) {
      rep("#9E9E9E", length(plot_obj$tip_cols))
    } else {
      plot_obj$tip_cols
    }
    internal_bg <- if (isTRUE(structure_only)) {
      rep("#9E9E9E", length(plot_obj$internal_cols))
    } else {
      plot_obj$internal_cols
    }
    tip_cex <- if (isTRUE(structure_only)) {
      rep(0.5, length(plot_obj$tip_cex))
    } else {
      plot_obj$tip_cex
    }
    internal_cex <- if (isTRUE(structure_only)) {
      rep(0.5, length(plot_obj$internal_cex))
    } else {
      plot_obj$internal_cex
    }
    ape::tiplabels(
      pch = 21,
      bg = tip_bg,
      col = "grey20",
      cex = tip_cex
    )
    ape::nodelabels(
      pch = 21,
      node = plot_obj$internal_node_ids,
      bg = internal_bg,
      col = "grey20",
      cex = internal_cex
    )
    if (!isTRUE(structure_only) && length(plot_obj$internal_ids) > 0L) {
      ape::nodelabels(
        text = plot_obj$internal_labels,
        node = plot_obj$internal_ids,
        frame = "none",
        cex = 0.7,
        font = 2,
        col = "black"
      )
    }
    if (!isTRUE(structure_only)) {
      legend(
        "topleft",
        inset = 0.01,
        bty = "n",
        title = "Representative sample",
        legend = c(
          paste0("Total counts: ", format(plot_obj$sample_total, big.mark = ",")),
          "Node fill = aggregated sample counts",
          "Labels = highest-count internal nodes"
        ),
        text.col = "grey20",
        cex = 0.85
      )
    }
  }

  grDevices::pdf(pdf_path, width = 14, height = 7)
  draw_once()
  grDevices::dev.off()

  grDevices::png(png_path, width = 14, height = 7, units = "in", res = 300)
  draw_once()
  grDevices::dev.off()
}

safe_scale <- function(x) {
  x_scaled <- scale(x)
  x_scaled[!is.finite(x_scaled)] <- 0
  x_scaled
}

clr_transform <- function(count_mat, pseudocount = 0.5) {
  x <- log(count_mat + pseudocount)
  sweep(x, 1, rowMeans(x), FUN = "-")
}

hellinger_transform <- function(count_mat) {
  row_total <- rowSums(count_mat)
  row_total[row_total <= 0] <- 1
  sqrt(count_mat / row_total)
}

split_logit_transform <- function(count_mat, tree, pseudocount = 0.5) {
  tree_info <- aggregate_tree_counts(count_mat, tree)
  agg <- as.matrix(tree_info$aggregated)
  sorted_nodes <- as.integer(tree_info$nodes)
  parent_nodes <- as.integer(tree_info$parent_nodes)
  node_to_col <- setNames(seq_along(sorted_nodes), sorted_nodes)
  edge_mat <- as.matrix(tree$edge)
  child_lookup <- split(edge_mat[, 2], edge_mat[, 1])

  split_mat <- matrix(NA_real_, nrow = nrow(agg), ncol = length(parent_nodes))
  colnames(split_mat) <- paste0("node_", parent_nodes)
  rownames(split_mat) <- rownames(count_mat)

  for (j in seq_along(parent_nodes)) {
    parent_node <- parent_nodes[[j]]
    child_nodes <- child_lookup[[as.character(parent_node)]]
    if (length(child_nodes) < 2L) {
      next
    }
    child_nodes <- sort(as.integer(child_nodes))
    parent_col <- node_to_col[[as.character(parent_node)]]
    left_col <- node_to_col[[as.character(child_nodes[[1L]])]]
    parent_count <- agg[, parent_col]
    left_count <- agg[, left_col]
    p_left <- (left_count + pseudocount) / (parent_count + 2 * pseudocount)
    split_mat[, j] <- qlogis(p_left)
  }

  split_mat[, colSums(is.finite(split_mat)) > 0, drop = FALSE]
}

compute_umap_embedding <- function(feature_mat) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop(
      "Package 'uwot' is required for UMAP embeddings. ",
      "Install it first with install.packages('uwot')."
    )
  }

  feature_mat <- as.matrix(feature_mat)
  storage.mode(feature_mat) <- "numeric"
  feature_mat[!is.finite(feature_mat)] <- 0
  feature_mat <- safe_scale(feature_mat)
  n_neighbors <- max(5L, min(30L, nrow(feature_mat) - 1L))
  embedding <- uwot::umap(
    feature_mat,
    n_neighbors = n_neighbors,
    min_dist = 0.2,
    metric = "euclidean",
    scale = FALSE,
    verbose = FALSE,
    ret_model = FALSE
  )
  colnames(embedding) <- c("Dim1", "Dim2")
  rownames(embedding) <- rownames(feature_mat)
  embedding
}

build_embedding_cluster_plot <- function(embedding, cluster, title, cluster_palette = NULL) {
  cluster_levels <- sort(unique(cluster))
  if (is.null(cluster_palette)) {
    cluster_palette <- setNames(
      scales::hue_pal(h = c(15, 375), c = 100, l = 65)(length(cluster_levels)),
      as.character(cluster_levels)
    )
  }

  plot_df <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    cluster = factor(cluster, levels = cluster_levels)
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, color = cluster)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.3) +
    ggplot2::scale_color_manual(values = cluster_palette) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Colored by CorTree cluster; n = ", nrow(plot_df)),
      x = "UMAP1",
      y = "UMAP2",
      color = "Cluster"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

build_embedding_covariate_plot <- function(embedding, cov_name, cov_df) {
  x_raw <- clean_covariate(cov_df[[cov_name]])
  x_num <- as_numeric_if_possible(x_raw)
  plot_df <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    stringsAsFactors = FALSE
  )

  if (!is.null(x_num)) {
    keep <- is.finite(x_num)
    plot_df <- cbind(plot_df[keep, , drop = FALSE], value = x_num[keep])
    return(
      ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, color = value)) +
        ggplot2::geom_point(alpha = 0.85, size = 1.3) +
        ggplot2::scale_color_gradient(low = "#D9E6F2", high = "#0B3C5D") +
        ggplot2::labs(
          title = cov_name,
          subtitle = "Colored by covariate value",
          x = "UMAP1",
          y = "UMAP2",
          color = cov_name
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 9),
          panel.grid.minor = ggplot2::element_blank()
        )
    )
  }

  x_chr <- as.character(x_raw)
  keep <- !is.na(x_chr)
  plot_df <- cbind(plot_df[keep, , drop = FALSE], value = x_chr[keep], stringsAsFactors = FALSE)
  value_count <- sort(table(plot_df$value), decreasing = TRUE)
  keep_levels <- names(value_count)[seq_len(min(10L, length(value_count)))]
  plot_df$value <- ifelse(plot_df$value %in% keep_levels, plot_df$value, "Other")
  plot_df$value <- factor(plot_df$value, levels = names(sort(table(plot_df$value), decreasing = TRUE)))

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Dim1, y = Dim2, color = value)) +
    ggplot2::geom_point(alpha = 0.85, size = 1.3) +
    ggplot2::labs(
      title = cov_name,
      subtitle = "Colored by covariate level",
      x = "UMAP1",
      y = "UMAP2",
      color = cov_name
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

build_numeric_cov_plot <- function(cov_name, cov_df, cluster_factor, method_label) {
  x_raw <- clean_covariate(cov_df[[cov_name]])
  x_num <- as_numeric_if_possible(x_raw)
  keep <- !is.na(x_num) & !is.na(cluster_factor)
  plot_df <- data.frame(
    cluster = cluster_factor[keep],
    value = x_num[keep]
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = cluster, y = value, fill = cluster)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.18, size = 0.8, color = "black") +
    ggplot2::labs(
      title = paste0("Top Numeric Covariate: ", cov_name),
      subtitle = method_label,
      x = "Final cluster",
      y = cov_name
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

build_categorical_cov_plot <- function(cov_name, cov_df, cluster_factor, method_label) {
  x_chr <- as.character(clean_covariate(cov_df[[cov_name]]))
  keep <- !is.na(x_chr) & !is.na(cluster_factor)
  plot_df <- data.frame(
    cluster = cluster_factor[keep],
    value = x_chr[keep],
    stringsAsFactors = FALSE
  )
  value_count <- sort(table(plot_df$value), decreasing = TRUE)
  keep_levels <- names(value_count)[seq_len(min(10L, length(value_count)))]
  plot_df$value <- ifelse(plot_df$value %in% keep_levels, plot_df$value, "Other")
  final_value_count <- sort(table(plot_df$value), decreasing = TRUE)
  final_levels <- names(final_value_count)
  plot_df$value <- factor(plot_df$value, levels = final_levels)
  value_labels <- setNames(
    paste0(final_levels, ":", as.integer(final_value_count[final_levels])),
    final_levels
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = cluster, fill = value)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::labs(
      title = cov_name,
      subtitle = method_label,
      x = "Final cluster",
      y = "Within-cluster proportion",
      fill = cov_name
    ) +
    ggplot2::scale_fill_discrete(labels = value_labels) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

build_covariate_summary_plot <- function(covariate_assoc, method_label) {
  if (is.null(covariate_assoc) || nrow(covariate_assoc) == 0L) {
    return(NULL)
  }

  top_assoc <- utils::head(covariate_assoc, 12L)
  top_assoc$covariate <- factor(top_assoc$covariate, levels = rev(top_assoc$covariate))
  ggplot2::ggplot(top_assoc, ggplot2::aes(x = effect_size, y = covariate, fill = var_type)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "Top covariates associated with final clustering",
      subtitle = method_label,
      x = "Effect size",
      y = NULL,
      fill = "Covariate type"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

build_pi_trace_plot <- function(pi_chain, iteration_index, cluster_factor, title, subtitle) {
  trace_df <- as.data.frame(t(pi_chain), check.names = FALSE)
  names(trace_df) <- paste0("Cluster ", seq.int(0L, ncol(trace_df) - 1L))
  trace_df$Iteration <- iteration_index
  trace_cols <- names(trace_df)[names(trace_df) != "Iteration"]
  trace_long <- stack(trace_df[trace_cols])
  names(trace_long) <- c("Value", "Cluster")
  trace_long$Iteration <- rep(trace_df$Iteration, times = length(trace_cols))
  trace_cluster_levels <- unique(as.character(trace_long$Cluster))
  cluster_size_table <- table(cluster_factor)
  trace_cluster_ids <- sub("^Cluster ", "", trace_cluster_levels)
  trace_cluster_sizes <- cluster_size_table[trace_cluster_ids]
  trace_cluster_sizes[is.na(trace_cluster_sizes)] <- 0L
  trace_labels <- setNames(
    paste0(trace_cluster_levels, ":", as.integer(trace_cluster_sizes)),
    trace_cluster_levels
  )
  trace_palette <- setNames(
    scales::hue_pal(h = c(15, 375), c = 100, l = 65)(length(trace_cluster_levels)),
    trace_cluster_levels
  )

  ggplot2::ggplot(trace_long, ggplot2::aes(x = Iteration, y = Value, color = Cluster)) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::scale_color_manual(values = trace_palette, labels = trace_labels) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "MCMC iteration index",
      y = expression(pi)
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_blank(),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )
}

build_dmm_weight_plot <- function(dmm_posterior, cluster_factor, mixture_weight = NULL) {
  mean_posterior <- colMeans(dmm_posterior, na.rm = TRUE)
  use_posterior_fallback <- is.null(mixture_weight) ||
    length(mixture_weight) != ncol(dmm_posterior)
  if (!use_posterior_fallback) {
    mixture_weight <- as.numeric(mixture_weight)
    use_posterior_fallback <- any(!is.finite(mixture_weight)) || sum(mixture_weight) <= 0
  }
  if (!use_posterior_fallback) {
    mixture_weight <- mixture_weight / sum(mixture_weight)
  }
  if (use_posterior_fallback) {
    mixture_weight <- mean_posterior / sum(mean_posterior)
  }
  cluster_ids <- seq_len(ncol(dmm_posterior)) - 1L
  weight_df <- data.frame(
    cluster = factor(cluster_ids, levels = cluster_ids),
    weight = as.numeric(mixture_weight),
    assigned_size = as.integer(table(factor(cluster_factor, levels = cluster_ids))),
    stringsAsFactors = FALSE
  )
  cluster_palette <- setNames(
    scales::hue_pal(h = c(15, 375), c = 100, l = 65)(length(cluster_ids)),
    as.character(cluster_ids)
  )
  cluster_labels <- setNames(
    paste0("cluster ", cluster_ids, ": ", weight_df$assigned_size),
    as.character(cluster_ids)
  )
  max_weight <- max(weight_df$weight, na.rm = TRUE)
  if (!is.finite(max_weight)) {
    max_weight <- 1
  }

  ggplot2::ggplot(weight_df, ggplot2::aes(x = cluster, y = weight, fill = cluster)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_fill_manual(
      values = cluster_palette,
      breaks = names(cluster_labels),
      labels = unname(cluster_labels),
      name = NULL
    ) +
    ggplot2::labs(
      title = "DMM cluster-weight diagnostic",
      subtitle = if (use_posterior_fallback) {
        "bars show posterior-average mixing proportions"
      } else {
        "bars show fitted mixing proportions"
      },
      x = "Cluster",
      y = "Mixing proportion",
      fill = NULL
    ) +
    ggplot2::coord_cartesian(
      ylim = c(0, min(1, max_weight + 0.08)),
      clip = "off"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
    )
}

compute_cluster_mean_ymax <- function(X, cluster, n_thresh = 10L) {
  cluster_chr <- as.character(round(cluster))
  cluster_counts <- table(cluster_chr)
  valid_clusters <- names(cluster_counts)[cluster_counts > n_thresh]

  if (length(valid_clusters) == 0L) {
    return(NA_real_)
  }

  profile_max <- vapply(
    valid_clusters,
    function(cluster_id) {
      max(colMeans(X[cluster_chr == cluster_id, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
    },
    numeric(1)
  )
  max(profile_max, na.rm = TRUE)
}

build_agp_cluster_mean_plot <- function(X, cluster, method_label, y_max = NULL, n_thresh = 10L) {
  cluster_mean_gg(
    X = X,
    chip_labels = cluster,
    y_max = y_max,
    n_thresh = n_thresh,
    title = method_label,
    x_label = "Microbiome taxa",
    y_label = "Count abundance"
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

build_agp_cluster_mean_overview <- function(
  X,
  cluster_assignments,
  method_titles = NULL,
  method_order = names(cluster_assignments),
  ncol = length(method_order),
  n_thresh = 10L
) {
  if (is.null(names(cluster_assignments)) || any(names(cluster_assignments) == "")) {
    stop("cluster_assignments must be a named list.")
  }

  method_order <- intersect(method_order, names(cluster_assignments))
  if (length(method_order) == 0L) {
    stop("No requested methods were found in cluster_assignments.")
  }

  if (is.null(method_titles)) {
    method_titles <- stats::setNames(method_order, method_order)
  }
  missing_titles <- setdiff(method_order, names(method_titles))
  if (length(missing_titles) > 0L) {
    method_titles[missing_titles] <- missing_titles
  }

  shared_y_max <- max(
    vapply(
      method_order,
      function(method_name) compute_cluster_mean_ymax(X, cluster_assignments[[method_name]], n_thresh = n_thresh),
      numeric(1)
    ),
    na.rm = TRUE
  )

  if (!is.finite(shared_y_max)) {
    stop("Could not compute a finite shared y-axis limit for the AGP cluster-mean overview.")
  }

  cluster_mean_plots <- lapply(method_order, function(method_name) {
    build_agp_cluster_mean_plot(
      X = X,
      cluster = cluster_assignments[[method_name]],
      method_label = unname(method_titles[[method_name]]),
      y_max = shared_y_max,
      n_thresh = n_thresh
    )
  })
  names(cluster_mean_plots) <- method_order

  patchwork::wrap_plots(cluster_mean_plots, ncol = ncol)
}
