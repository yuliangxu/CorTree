rm(list = ls())

project_dir <- normalizePath(".")
source(file.path(project_dir, "R", "AGP_help.R"))

library(ggplot2)
library(tidyr)
library(patchwork)

if (!requireNamespace("ape", quietly = TRUE)) {
  stop(
    "AGP2_analysis.R needs package 'ape' to draw the phylogenetic tree figure. ",
    "Please run install.packages('ape') and then rerun the script."
  )
}

out_dir <- "/cwork/yx306/CorTree"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bundle_path_candidates <- file.path(out_dir, "AGP1_phylotree_fit.rds")
bundle_path_candidates <- unique(bundle_path_candidates[file.exists(bundle_path_candidates)])
if (length(bundle_path_candidates) == 0L) {
  stop("Could not find AGP1 result bundle: AGP1_phylotree_fit.rds")
}
bundle_path <- bundle_path_candidates[[1L]]

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
ag_cov <- ag_obj$ag_fecal
cov_clean_info <- clean_covariate_dataframe(ag_cov)
ag_cov_clean <- cov_clean_info$covariates
cov_normalization_report <- cov_clean_info$normalization_report
sampler_args <- bundle$input$sampler_args
fit <- bundle$fit
tree <- ag_obj$tree

table(ag_cov_clean$COUNTRY_RESIDENCE)
summary(ag_cov_clean$COUNTRY_RESIDENCE)

if (is.null(keep_sample) || length(keep_sample) != nrow(count_data)) {
  keep_sample <- rownames(count_data)
}

Z_post <- fit$mcmc$Z
cluster_hat <- apply(Z_post, 1, function(z) {
  as.integer(names(which.max(table(z))))
})

covariate_assoc <- analyze_covariate_association(
  cluster = cluster_hat,
  covariates = ag_cov_clean,
  sample_ids = keep_sample
)

top_covariate <- if (nrow(covariate_assoc) > 0L) covariate_assoc[1, , drop = FALSE] else NULL
top_categorical <- if (nrow(covariate_assoc) > 0L) {
  covariate_assoc[covariate_assoc$var_type == "categorical", , drop = FALSE]
} else {
  data.frame()
}
top_categorical_1 <- if (nrow(top_categorical) > 0L) top_categorical[1, , drop = FALSE] else NULL
top_categorical_2 <- if (nrow(top_categorical) > 1L) top_categorical[2, , drop = FALSE] else NULL

cov_df <- as.data.frame(ag_cov_clean, stringsAsFactors = FALSE)
row_idx <- match(keep_sample, cov_df$SampleID)
cov_df <- cov_df[row_idx, , drop = FALSE]
cluster_factor <- factor(cluster_hat, levels = sort(unique(cluster_hat)))

build_numeric_cov_plot <- function(cov_name, cov_df, cluster_factor) {
  x_raw <- clean_covariate(cov_df[[cov_name]])
  x_num <- as_numeric_if_possible(x_raw)
  keep <- !is.na(x_num) & !is.na(cluster_factor)
  plot_df <- data.frame(
    cluster = cluster_factor[keep],
    value = x_num[keep]
  )
  ggplot(plot_df, aes(x = cluster, y = value, fill = cluster)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.18, size = 0.8, color = "black") +
    labs(
      title = paste0("Top Numeric Covariate: ", cov_name),
      x = "Final cluster",
      y = cov_name
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

build_categorical_cov_plot <- function(cov_name, cov_df, cluster_factor) {
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

  ggplot(plot_df, aes(x = cluster, fill = value)) +
    geom_bar(position = "fill") +
    labs(
      title = paste0(cov_name),
      x = "Final cluster",
      y = "Within-cluster proportion",
      fill = cov_name
    ) +
    scale_fill_discrete(labels = value_labels) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

pi_trace_info <- get_pi_trace_plot_data(fit = fit, sampler_args = sampler_args)
pi_chain <- pi_trace_info$pi_chain
pi_iteration_index <- pi_trace_info$iteration_index

trace_df <- as.data.frame(t(pi_chain), check.names = FALSE)
names(trace_df) <- paste0("Cluster ", seq.int(0L, ncol(trace_df) - 1L))
trace_df$Iteration <- pi_iteration_index
trace_long <- tidyr::pivot_longer(
  trace_df,
  cols = starts_with("Cluster "),
  names_to = "Cluster",
  values_to = "Value"
)
trace_cluster_levels <- unique(trace_long$Cluster)
cluster_size_table <- table(cluster_factor)
trace_cluster_ids <- sub("^Cluster ", "", trace_cluster_levels)
trace_cluster_sizes <- cluster_size_table[trace_cluster_ids]
trace_cluster_sizes[is.na(trace_cluster_sizes)] <- 0L
trace_labels <- setNames(
  paste0(
    trace_cluster_levels,
    ":",
    as.integer(trace_cluster_sizes)
  ),
  trace_cluster_levels
)
trace_palette <- setNames(
  scales::hue_pal(h = c(15, 375), c = 100, l = 65)(length(trace_cluster_levels)),
  trace_cluster_levels
)

pi_trace_plot <- ggplot(trace_long, aes(x = Iteration, y = Value, color = Cluster)) +
  geom_line(linewidth = 0.45) +
  scale_color_manual(values = trace_palette, labels = trace_labels) +
  labs(
    title = "AGP mixing-weight trace plot",
    subtitle = if (identical(pi_trace_info$trace_source, "full")) {
      paste0("Full pi trace from iteration ", min(pi_iteration_index), " to ", max(pi_iteration_index))
    } else if (identical(pi_trace_info$trace_source, "post_burnin_only")) {
      paste0("Older fit bundle stores post-burnin pi only; plotted iterations are labeled ", min(pi_iteration_index), "-", max(pi_iteration_index))
    } else {
      paste0("Using saved pi trace labeled ", min(pi_iteration_index), "-", max(pi_iteration_index))
    },
    x = "MCMC iteration index",
    y = expression(pi)
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(size = 9),
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

cov_summary_plot <- NULL
if (nrow(covariate_assoc) > 0L) {
  top_assoc <- utils::head(covariate_assoc, 12L)
  top_assoc$covariate <- factor(top_assoc$covariate, levels = rev(top_assoc$covariate))
  cov_summary_plot <- ggplot(top_assoc, aes(x = effect_size, y = covariate, fill = var_type)) +
    geom_col() +
    labs(
      title = "Top covariates associated with final clustering",
      x = "Effect size",
      y = NULL,
      fill = "Covariate type"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

demo_plots <- list(pi_trace_plot)
if (!is.null(top_categorical_1)) {
  demo_plots[[length(demo_plots) + 1L]] <- build_categorical_cov_plot(top_categorical_1$covariate[[1L]], cov_df, cluster_factor)
}
if (!is.null(top_categorical_2)) {
  demo_plots[[length(demo_plots) + 1L]] <- build_categorical_cov_plot(top_categorical_2$covariate[[1L]], cov_df, cluster_factor)
}
if (!is.null(cov_summary_plot)) {
  demo_plots[[length(demo_plots) + 1L]] <- cov_summary_plot
}
combined_plot <- wrap_plots(demo_plots, ncol = 2)

phylo_plot_obj <- NULL
if (inherits(tree, "phylo")) {
  phylo_plot_obj <- build_phylo_node_plot(
    tree = tree,
    count_data = count_data
  )
}

trace_plot_pdf <- file.path(out_dir, "AGP2_pi_trace_full_ggplot.pdf")
trace_plot_png <- file.path(out_dir, "AGP2_pi_trace_full_ggplot.png")
cov_plot_pdf <- file.path(out_dir, "AGP2_cluster_covariate_diagnostics.pdf")
cov_plot_png <- file.path(out_dir, "AGP2_cluster_covariate_diagnostics.png")
phylo_plot_pdf <- file.path(out_dir, "AGP2_phylogenetic_tree_sample_counts.pdf")
phylo_plot_png <- file.path(out_dir, "AGP2_phylogenetic_tree_sample_counts.png")
cov_assoc_out <- file.path(out_dir, "AGP2_cluster_covariate_association.csv")
cluster_out <- file.path(out_dir, "AGP2_cluster_assignments.csv")
cov_clean_out <- file.path(out_dir, "AGP2_covariates_cleaned.csv")
cov_clean_rds_out <- file.path(out_dir, "AGP2_covariates_cleaned.rds")
cov_clean_report_out <- file.path(out_dir, "AGP2_covariate_numeric_label_normalization.csv")
summary_txt_out <- file.path(out_dir, "AGP2_analysis_summary.txt")
plan_txt_out <- file.path(out_dir, "AGP2_analysis_plan.txt")
analysis_rds_out <- file.path(out_dir, "AGP2_analysis_bundle.rds")

ggplot2::ggsave(filename = trace_plot_pdf, plot = pi_trace_plot, width = 10, height = 4.5)
ggplot2::ggsave(filename = trace_plot_png, plot = pi_trace_plot, width = 10, height = 4.5, dpi = 300)
ggplot2::ggsave(filename = cov_plot_pdf, plot = combined_plot, width = 14, height = 9)
ggplot2::ggsave(filename = cov_plot_png, plot = combined_plot, width = 14, height = 9, dpi = 300)
if (!is.null(phylo_plot_obj)) {
  save_phylo_node_plot(
    plot_obj = phylo_plot_obj,
    pdf_path = phylo_plot_pdf,
    png_path = phylo_plot_png,
    structure_only = TRUE
  )
}

utils::write.csv(covariate_assoc, cov_assoc_out, row.names = FALSE)
utils::write.csv(
  data.frame(
    SampleID = keep_sample,
    cluster_hat = cluster_hat,
    stringsAsFactors = FALSE
  ),
  cluster_out,
  row.names = FALSE
)
utils::write.csv(cov_df, cov_clean_out, row.names = FALSE, na = "")
saveRDS(cov_df, cov_clean_rds_out)
utils::write.csv(cov_normalization_report, cov_clean_report_out, row.names = FALSE)

plan_lines <- c(
  "Plan to demonstrate how observed covariates relate to the unsupervised clustering result",
  "",
  "1. Start with an unsupervised summary.",
  "Use the final clustering assignment from the posterior MCMC samples and report cluster sizes so the covariate analysis is interpreted relative to the discovered groups rather than imposed labels.",
  "",
  "2. Rank every observed covariate with a type-aware association test.",
  "For numeric covariates use Kruskal-Wallis and report epsilon-squared; for categorical covariates use chi-squared and report Cramer's V. This gives one comparable screening table across all columns in ag_cov.",
  "",
  "3. Focus interpretation on effect size first, p-value second.",
  "The most informative covariates are the ones with the largest effect sizes and acceptable adjusted p-values, not just the smallest raw p-values in a large table.",
  "",
  "4. Visualize the top signals in cluster-conditional form.",
  "For top numeric covariates, show boxplots or violin plots by cluster. For top categorical covariates, show within-cluster composition bars so we can see whether a covariate is separating one cluster strongly or shifting all clusters gradually.",
  "",
  "5. Distinguish descriptive association from causal interpretation.",
  "Because the clustering is unsupervised, these plots show which metadata align with the latent structure; they do not prove that the covariates caused the clustering or that the model used those covariates directly.",
  "",
  "6. Check robustness of the interpretation.",
  "If we want stronger evidence, the next step is to compare the same covariate rankings across multiple random seeds or initialization modes, or to use posterior uncertainty in Z rather than only the final hard clustering.",
  "",
  "7. Important limitation for the current pi trace plot.",
  paste0(
    "The current AGP fit stores only ",
    ncol(pi_chain),
    " post-burnin pi draws, so the trace plot covers saved iterations ",
    min(pi_iteration_index),
    "-",
    max(pi_iteration_index),
    " rather than the full warmup-plus-sampling chain."
  )
)

summary_lines <- c(
  "AGP2 analysis summary",
  "",
  paste("Bundle path:", bundle_path),
  paste("Data path:", data_path),
  paste("Samples:", nrow(count_data)),
  paste("Taxa:", ncol(count_data)),
  paste("Saved pi draws available:", ncol(pi_chain)),
  if (!is.null(phylo_plot_obj)) {
    paste(
      "Tree figure sample:",
      phylo_plot_obj$sample_id,
      "(total counts =",
      format(phylo_plot_obj$sample_total, big.mark = ","),
      ")"
    )
  } else {
    "Tree figure sample: tree object is not of class phylo, so no tree plot was generated."
  },
  paste("Cluster sizes:", paste(names(table(cluster_hat)), table(cluster_hat), sep = "=", collapse = ", ")),
  paste("Cleaned covariate file:", cov_clean_out),
  paste("Normalized numeric-label issues found:", nrow(cov_normalization_report)),
  if (nrow(cov_normalization_report) > 0L) {
    paste(
      "Affected covariates:",
      paste(unique(cov_normalization_report$covariate), collapse = ", ")
    )
  } else {
    "Affected covariates: none"
  },
  "",
  "Top covariates by association effect size:",
  capture.output(print(utils::head(covariate_assoc, 10L), row.names = FALSE)),
  "",
  "Most informative covariate:",
  if (!is.null(top_covariate)) {
    capture.output(print(top_covariate, row.names = FALSE))
  } else {
    "No eligible covariate association was found."
  }
)

writeLines(plan_lines, con = plan_txt_out)
writeLines(summary_lines, con = summary_txt_out)

saveRDS(
  list(
    bundle_path = bundle_path,
    data_path = data_path,
    cluster_hat = cluster_hat,
    covariates_clean = cov_df,
    covariate_normalization_report = cov_normalization_report,
    covariate_association = covariate_assoc,
    top_covariate = top_covariate,
    phylo_tree_plot = if (!is.null(phylo_plot_obj)) {
      list(
        sample_id = phylo_plot_obj$sample_id,
        sample_total = phylo_plot_obj$sample_total,
        plot_pdf = phylo_plot_pdf,
        plot_png = phylo_plot_png
      )
    } else {
      NULL
    },
    pi_iteration_index = pi_iteration_index,
    plan = plan_lines
  ),
  analysis_rds_out
)

cat("AGP2 analysis complete\n")
cat("Saved pi trace plot:", trace_plot_pdf, "\n")
cat("Saved pi trace plot:", trace_plot_png, "\n")
cat("Saved covariate diagnostics:", cov_plot_pdf, "\n")
cat("Saved covariate diagnostics:", cov_plot_png, "\n")
if (!is.null(phylo_plot_obj)) {
  cat("Saved phylogenetic tree figure:", phylo_plot_pdf, "\n")
  cat("Saved phylogenetic tree figure:", phylo_plot_png, "\n")
}
cat("Saved covariate association table:", cov_assoc_out, "\n")
cat("Saved cluster assignments:", cluster_out, "\n")
cat("Saved cleaned covariates:", cov_clean_out, "\n")
cat("Saved cleaned covariates RDS:", cov_clean_rds_out, "\n")
cat("Saved covariate normalization report:", cov_clean_report_out, "\n")
cat("Saved text summary:", summary_txt_out, "\n")
cat("Saved analysis plan:", plan_txt_out, "\n")
cat("Saved analysis bundle:", analysis_rds_out, "\n")
if (!is.null(top_covariate)) {
  cat("Most informative covariate:\n")
  print(top_covariate)
}
