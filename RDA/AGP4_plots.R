rm(list = ls())

project_dir <- normalizePath(".")
devtools::load_all()
source(file.path(project_dir, "R", "AGP_help.R"))

library(ggplot2)
library(patchwork)

if (!requireNamespace("uwot", quietly = TRUE)) {
  stop(
    "AGP4_plots.R needs package 'uwot'. ",
    "Please install it before running the script."
  )
}

out_dir <- "/cwork/yx306/CorTree"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_gc <- function(...) {
  invisible(gc(verbose = FALSE, ...))
}

bundle_path_candidates <- unique(c(
  file.path(out_dir, "AGP1_phylotree_fit.rds")
))
bundle_path_candidates <- bundle_path_candidates[file.exists(bundle_path_candidates)]
if (length(bundle_path_candidates) == 0L) {
  stop("Could not find AGP1 result bundle: AGP1_phylotree_fit.rds")
}
bundle_path <- bundle_path_candidates[[1L]]

compare_bundle_candidates <- unique(c(
  file.path(out_dir, "AGP3_compare_bundle.rds")
))
compare_bundle_candidates <- compare_bundle_candidates[file.exists(compare_bundle_candidates)]
compare_bundle_path <- if (length(compare_bundle_candidates) > 0L) {
  compare_bundle_candidates[[1L]]
} else {
  NULL
}

bundle <- readRDS(bundle_path)
if (is.null(bundle$fit) || is.null(bundle$input)) {
  stop("Expected AGP1_phylotree_fit.rds to contain the AGP1 result bundle with $fit and $input.")
}

data_path <- bundle$input$data_path
if (!file.exists(data_path)) {
  stop("Data path recorded in the AGP1 bundle does not exist: ", data_path)
}

ag_obj <- readRDS(data_path)
count_data <- t(as.matrix(ag_obj$otu_top))
keep_sample <- bundle$input$keep_sample
if (is.null(keep_sample) || length(keep_sample) != nrow(count_data)) {
  keep_sample <- rownames(count_data)
}

fit <- bundle$fit
Z_post <- fit$mcmc$Z
cluster_hat <- apply(Z_post, 1, function(z) {
  as.integer(names(which.max(table(z))))
})
cluster_levels <- sort(unique(cluster_hat))
cluster_palette <- setNames(
  scales::hue_pal(h = c(15, 375), c = 100, l = 65)(length(cluster_levels)),
  as.character(cluster_levels)
)

ag_cov <- ag_obj$ag_fecal
cov_clean_info <- clean_covariate_dataframe(ag_cov)
ag_cov_clean <- cov_clean_info$covariates
covariate_assoc <- analyze_covariate_association(
  cluster = cluster_hat,
  covariates = ag_cov_clean,
  sample_ids = keep_sample
)

cov_df <- as.data.frame(ag_cov_clean, stringsAsFactors = FALSE)
row_idx <- match(keep_sample, cov_df$SampleID)
cov_df <- cov_df[row_idx, , drop = FALSE]
ag_tree <- ag_obj$tree
rm(ag_obj, ag_cov, cov_clean_info, ag_cov_clean, row_idx)
run_gc()

embedding_title <- c(
  log1p = "UMAP of log1p(count_data)",
  clr = "UMAP of CLR(count_data + 0.5)",
  hellinger = "UMAP of Hellinger(count_data)",
  split_logit = "UMAP of empirical logit split probabilities"
)

embedding_builders <- list(
  log1p = function() log1p(count_data),
  clr = function() clr_transform(count_data, pseudocount = 0.5),
  hellinger = function() hellinger_transform(count_data),
  split_logit = function() split_logit_transform(count_data, ag_tree, pseudocount = 0.5)
)

embedding_list <- setNames(vector("list", length(embedding_builders)), names(embedding_builders))
for (name in names(embedding_builders)) {
  feature_matrix <- embedding_builders[[name]]()
  embedding_list[[name]] <- compute_umap_embedding(feature_matrix)
  rm(feature_matrix)
  run_gc()
}
rm(embedding_builders, ag_tree)
run_gc()

cluster_plot_list <- lapply(names(embedding_list), function(name) {
  build_embedding_cluster_plot(
    embedding = embedding_list[[name]],
    cluster = cluster_hat,
    title = embedding_title[[name]],
    cluster_palette = cluster_palette
  )
})
names(cluster_plot_list) <- names(embedding_list)

cluster_overview <- wrap_plots(cluster_plot_list, ncol = 2)

top_covariates <- if (nrow(covariate_assoc) > 0L) {
  utils::head(covariate_assoc$covariate, 3L)
} else {
  character()
}

covariate_plot_files <- list()
if (length(top_covariates) > 0L) {
  for (cov_name in top_covariates) {
    cov_plot_list <- lapply(names(embedding_list), function(name) {
      build_embedding_covariate_plot(
        embedding = embedding_list[[name]],
        cov_name = cov_name,
        cov_df = cov_df
      ) +
        ggtitle(paste0(embedding_title[[name]], "\n", cov_name))
    })
    cov_combined <- wrap_plots(cov_plot_list, ncol = 2)
    cov_file_stub <- gsub("[^A-Za-z0-9]+", "_", cov_name)
    cov_png <- file.path(out_dir, paste0("AGP4_covariate_", cov_file_stub, "_overview.png"))
    ggplot2::ggsave(cov_png, cov_combined, width = 14, height = 10, dpi = 300)
    covariate_plot_files[[cov_name]] <- cov_png
    rm(cov_plot_list, cov_combined)
    run_gc()
  }
}

cluster_overview_png <- file.path(out_dir, "AGP4_umap_cluster_overview.png")
cluster_overview_pdf <- file.path(out_dir, "AGP4_umap_cluster_overview.pdf")
method_comparison_png <- file.path(out_dir, "AGP4_real_data_method_comparison.png")
method_comparison_pdf <- file.path(out_dir, "AGP4_real_data_method_comparison.pdf")
cluster_mean_comparison_png <- file.path(out_dir, "AGP4_real_data_cluster_mean_comparison.png")
cluster_mean_comparison_pdf <- file.path(out_dir, "AGP4_real_data_cluster_mean_comparison.pdf")
summary_txt_out <- file.path(out_dir, "AGP4_plot_summary.txt")
plot_bundle_out <- file.path(out_dir, "AGP4_plot_bundle.rds")
embedding_csv_out <- file.path(out_dir, "AGP4_umap_embeddings.csv")

ggplot2::ggsave(cluster_overview_png, cluster_overview, width = 14, height = 10, dpi = 300)
ggplot2::ggsave(cluster_overview_pdf, cluster_overview, width = 14, height = 10)

embedding_df <- do.call(
  rbind,
  lapply(names(embedding_list), function(name) {
    embedding <- embedding_list[[name]]
    data.frame(
      SampleID = rownames(embedding),
      embedding = name,
      Dim1 = embedding[, 1],
      Dim2 = embedding[, 2],
      cluster_hat = cluster_hat,
      stringsAsFactors = FALSE
    )
  })
)
utils::write.csv(embedding_df, embedding_csv_out, row.names = FALSE)
rm(embedding_df)
run_gc()


# comparison with AGP3 methods if comparison bundle is available ===================================================
devtools::load_all()
source(file.path(project_dir, "R", "AGP_help.R"))
comparison_outputs <- list()
comparison_summary_lines <- "Real-data method comparison: AGP3 comparison bundle not found; skipped comparison diagnostics."
if (!is.null(compare_bundle_path)) {
  compare_bundle <- readRDS(compare_bundle_path)
  compare_assignments <- compare_bundle$cluster_assignments
  if (!is.null(compare_assignments) && nrow(compare_assignments) == length(keep_sample)) {
    compare_idx <- match(keep_sample, compare_assignments$SampleID)
    if (anyNA(compare_idx)) {
      stop("AGP3 comparison bundle sample IDs do not align with AGP4 samples.")
    }
    compare_assignments <- compare_assignments[compare_idx, , drop = FALSE]
  }

  compare_umap <- compare_bundle$umap
  if (is.null(compare_umap)) {
    compare_umap <- embedding_list$log1p
  }
  if (
    is.null(compare_assignments) ||
    !all(c("CorTree", "DMM", "IndTree") %in% names(compare_assignments))
  ) {
    comparison_summary_lines <- c(
      paste("Real-data comparison bundle:", compare_bundle_path),
      "Comparison diagnostics skipped: cluster assignments are missing required method columns."
    )
  } else {
    method_clusters <- list(
      CorTree = compare_assignments$CorTree,
      DMM = compare_assignments$DMM,
      IndTree = compare_assignments$IndTree
    )
    method_cluster_mean_plot <- build_agp_cluster_mean_overview(
      X = count_data,
      cluster_assignments = method_clusters,
      method_titles = c(
        CorTree = "CorTree",
        DMM = "DMM",
        IndTree = "IndTree"
      ),
      method_order = c("CorTree", "DMM", "IndTree"),
      ncol = 3
    )
    ggplot2::ggsave(
      cluster_mean_comparison_png,
      method_cluster_mean_plot,
      width = 16,
      height = 5.5,
      dpi = 300
    )
    ggplot2::ggsave(
      cluster_mean_comparison_pdf,
      method_cluster_mean_plot,
      width = 16,
      height = 5.5
    )
    method_palette <- lapply(method_clusters, function(cluster_vec) {
      levels_chr <- as.character(sort(unique(cluster_vec)))
      setNames(
        scales::hue_pal(h = c(15, 375), c = 100, l = 65)(length(levels_chr)),
        levels_chr
      )
    })
    ari_lookup <- if (!is.null(compare_bundle$method_summary)) {
      setNames(compare_bundle$method_summary$ari_vs_cortree, compare_bundle$method_summary$method)
    } else {
      c(CorTree = 1, DMM = NA_real_, IndTree = NA_real_)
    }
    comparison_cluster_plots <- lapply(names(method_clusters), function(method_name) {
      cluster_vec <- method_clusters[[method_name]]
      subtitle_txt <- if (identical(method_name, "CorTree")) {
        paste0("Reference clustering; n = ", length(cluster_vec))
      } else {
        paste0("ARI vs CorTree = ", sprintf("%.3f", ari_lookup[[method_name]]), "; n = ", length(cluster_vec))
      }
      build_embedding_cluster_plot(
        embedding = compare_umap,
        cluster = cluster_vec,
        title = paste0("AGP clustering on shared UMAP: ", method_name),
        cluster_palette = method_palette[[method_name]]
      ) +
        labs(subtitle = subtitle_txt)
    })
    names(comparison_cluster_plots) <- names(method_clusters)
    method_comparison_plot <- wrap_plots(comparison_cluster_plots, ncol = 3)
    ggplot2::ggsave(method_comparison_png, method_comparison_plot, width = 16, height = 5.5, dpi = 300)
    ggplot2::ggsave(method_comparison_pdf, method_comparison_plot, width = 16, height = 5.5)
    rm(comparison_cluster_plots, method_comparison_plot)
    run_gc()

    cov_assoc_by_method <- compare_bundle$covariate_association
    method_diagnostic_files <- list()
    comparison_summary_lines <- c(
      paste("Real-data comparison bundle:", compare_bundle_path),
      paste("Comparison methods:", paste(names(method_clusters), collapse = ", ")),
      paste("Cluster-mean overview PNG:", cluster_mean_comparison_png),
      paste("Cluster-mean overview PDF:", cluster_mean_comparison_pdf),
      paste(
        "ARI vs CorTree:",
        paste(
          paste0(names(ari_lookup), "=", sprintf("%.3f", ari_lookup)),
          collapse = ", "
        )
      )
    )

    cortree_pi_trace_info <- get_pi_trace_plot_data(
      fit = fit,
      sampler_args = bundle$input$sampler_args
    )
    indtree_pi_trace_info <- NULL
    if (!is.null(compare_bundle$fit$indtree)) {
      indtree_pi_trace_info <- get_pi_trace_plot_data(
        fit = compare_bundle$fit$indtree,
        sampler_args = compare_bundle$input$sampler_args
      )
    } else if (!is.null(compare_bundle$trace$indtree_pi_postburnin)) {
      indtree_pi_trace <- as.matrix(compare_bundle$trace$indtree_pi_postburnin)
      storage.mode(indtree_pi_trace) <- "numeric"
      indtree_burnin <- if (!is.null(compare_bundle$input$sampler_args$burnin)) {
        as.integer(compare_bundle$input$sampler_args$burnin)
      } else {
        0L
      }
      indtree_pi_trace_info <- list(
        pi_chain = indtree_pi_trace,
        iteration_index = seq.int(indtree_burnin, indtree_burnin + ncol(indtree_pi_trace) - 1L),
        trace_source = "post_burnin_only"
      )
    }
    dmm_posterior <- compare_bundle$posterior$dmm
    dmm_mixture_weight <- compare_bundle$metadata$dmm_mixture_weight

    for (method_name in names(method_clusters)) {
      cluster_vec <- as.integer(method_clusters[[method_name]])
      cluster_factor_method <- factor(cluster_vec, levels = sort(unique(cluster_vec)))
      assoc_key <- tolower(method_name)
      method_assoc <- cov_assoc_by_method[[assoc_key]]

      main_plot <- NULL
      if (identical(method_name, "CorTree")) {
        main_plot <- build_pi_trace_plot(
          pi_chain = cortree_pi_trace_info$pi_chain,
          iteration_index = cortree_pi_trace_info$iteration_index,
          cluster_factor = cluster_factor_method,
          title = "CorTree mixing-weight trace plot",
          subtitle = "Same diagnostic style as AGP2_analysis.R"
        )
      } else if (identical(method_name, "IndTree") && !is.null(indtree_pi_trace_info)) {
        main_plot <- build_pi_trace_plot(
          pi_chain = indtree_pi_trace_info$pi_chain,
          iteration_index = indtree_pi_trace_info$iteration_index,
          cluster_factor = cluster_factor_method,
          title = "IndTree mixing-weight trace plot",
          subtitle = "Same diagnostic style as AGP2_analysis.R"
        )
      } else if (identical(method_name, "DMM") && !is.null(dmm_posterior)) {
        main_plot <- build_dmm_weight_plot(
          dmm_posterior = dmm_posterior,
          cluster_factor = cluster_vec,
          mixture_weight = dmm_mixture_weight
        )
      }

      method_plots <- Filter(Negate(is.null), list(main_plot))
      top_categorical <- if (!is.null(method_assoc) && nrow(method_assoc) > 0L) {
        method_assoc[method_assoc$var_type == "categorical", , drop = FALSE]
      } else {
        data.frame()
      }
      top_numeric <- if (!is.null(method_assoc) && nrow(method_assoc) > 0L) {
        method_assoc[method_assoc$var_type == "numeric", , drop = FALSE]
      } else {
        data.frame()
      }

      if (nrow(top_categorical) > 0L) {
        method_plots[[length(method_plots) + 1L]] <- build_categorical_cov_plot(
          top_categorical$covariate[[1L]],
          cov_df,
          cluster_factor_method,
          method_name
        )
      }
      if (nrow(top_categorical) > 1L) {
        method_plots[[length(method_plots) + 1L]] <- build_categorical_cov_plot(
          top_categorical$covariate[[2L]],
          cov_df,
          cluster_factor_method,
          method_name
        )
      } else if (nrow(top_numeric) > 0L) {
        method_plots[[length(method_plots) + 1L]] <- build_numeric_cov_plot(
          top_numeric$covariate[[1L]],
          cov_df,
          cluster_factor_method,
          method_name
        )
      }

      cov_summary_plot <- build_covariate_summary_plot(method_assoc, method_name)
      if (!is.null(cov_summary_plot)) {
        method_plots[[length(method_plots) + 1L]] <- cov_summary_plot
      }

      if (length(method_plots) == 0L) {
        next
      }

      method_diagnostic_plot <- wrap_plots(method_plots, ncol = 2)
      method_file_stub <- tolower(method_name)
      method_png <- file.path(out_dir, paste0("AGP4_", method_file_stub, "_diagnostic_overview.png"))
      method_pdf <- file.path(out_dir, paste0("AGP4_", method_file_stub, "_diagnostic_overview.pdf"))
      ggplot2::ggsave(method_png, method_diagnostic_plot, width = 14, height = 9, dpi = 300)
      ggplot2::ggsave(method_pdf, method_diagnostic_plot, width = 14, height = 9)

      method_diagnostic_files[[method_name]] <- c(png = method_png, pdf = method_pdf)
      comparison_summary_lines <- c(
        comparison_summary_lines,
        paste0(
          method_name,
          " top covariates: ",
          if (!is.null(method_assoc) && nrow(method_assoc) > 0L) {
            paste(utils::head(method_assoc$covariate, 3L), collapse = ", ")
          } else {
            "none detected"
          }
        )
      )
      rm(
        cluster_vec,
        cluster_factor_method,
        method_assoc,
        main_plot,
        method_plots,
        top_categorical,
        top_numeric,
        cov_summary_plot,
        method_diagnostic_plot
      )
      run_gc()
    }

    comparison_outputs <- list(
      compare_bundle_path = compare_bundle_path,
      cluster_mean_comparison_png = cluster_mean_comparison_png,
      cluster_mean_comparison_pdf = cluster_mean_comparison_pdf,
      method_comparison_png = method_comparison_png,
      method_comparison_pdf = method_comparison_pdf,
      method_diagnostic_files = method_diagnostic_files
    )
    rm(
      compare_bundle,
      compare_assignments,
      compare_umap,
      method_cluster_mean_plot,
      method_clusters,
      method_palette,
      ari_lookup,
      cov_assoc_by_method,
      method_diagnostic_files,
      cortree_pi_trace_info,
      indtree_pi_trace_info,
      dmm_posterior,
      dmm_mixture_weight
    )
    run_gc()
  }
}

summary_lines <- c(
  "AGP4 candidate embedding plots for CorTree result",
  "",
  paste("Bundle path:", bundle_path),
  paste("Data path:", data_path),
  paste("Samples:", nrow(count_data)),
  paste("Taxa:", ncol(count_data)),
  paste("Cluster sizes:", paste(names(table(cluster_hat)), table(cluster_hat), sep = "=", collapse = ", ")),
  "",
  "Candidate UMAP inputs:",
  "1. log1p(count_data)",
  "2. CLR(count_data + 0.5)",
  "3. Hellinger(count_data)",
  "4. empirical logit split probabilities from the phylogenetic tree",
  "",
  "Interpretation note:",
  "The split-logit embedding uses empirical node splitting proportions from the tree, with the first child chosen by node order. This is useful for exploratory diagnostics, but the sign of a split can depend on child ordering.",
  "",
  comparison_summary_lines,
  "",
  if (length(top_covariates) > 0L) {
    paste("Top covariates plotted:", paste(top_covariates, collapse = ", "))
  } else {
    "Top covariates plotted: none detected from the current association table."
  }
)
writeLines(summary_lines, summary_txt_out)

saveRDS(
  list(
    input = list(
      bundle_path = bundle_path,
      data_path = data_path
    ),
    cluster_hat = cluster_hat,
    covariate_association = covariate_assoc,
    top_covariates = top_covariates,
    embeddings = embedding_list,
    feature_metadata = list(
      transformations = names(embedding_title),
      split_logit_pseudocount = 0.5,
      clr_pseudocount = 0.5
    ),
    output = list(
      cluster_overview_png = cluster_overview_png,
      cluster_overview_pdf = cluster_overview_pdf,
      covariate_plot_files = covariate_plot_files,
      summary_txt = summary_txt_out,
      embedding_csv = embedding_csv_out
    ),
    comparison = comparison_outputs
  ),
  plot_bundle_out
)

cat("AGP4 plot generation complete\n")
cat("Saved cluster overview PNG:", cluster_overview_png, "\n")
cat("Saved cluster overview PDF:", cluster_overview_pdf, "\n")
cat("Saved embedding CSV:", embedding_csv_out, "\n")
cat("Saved plot bundle:", plot_bundle_out, "\n")
cat("Saved summary:", summary_txt_out, "\n")
if (length(comparison_outputs) > 0L) {
  cat("Saved cluster mean comparison PNG:", comparison_outputs$cluster_mean_comparison_png, "\n")
  cat("Saved cluster mean comparison PDF:", comparison_outputs$cluster_mean_comparison_pdf, "\n")
  cat("Saved method comparison PNG:", comparison_outputs$method_comparison_png, "\n")
  cat("Saved method comparison PDF:", comparison_outputs$method_comparison_pdf, "\n")
  cat("Saved method diagnostic plots:\n")
  print(comparison_outputs$method_diagnostic_files)
}
if (length(covariate_plot_files) > 0L) {
  cat("Saved covariate overview PNGs:\n")
  print(unlist(covariate_plot_files, use.names = TRUE))
}
