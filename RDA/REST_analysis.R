rm(list = ls())

project_dir <- normalizePath(".")
source(file.path(project_dir, "R", "DCC_help.R"))
source(file.path(project_dir, "R", "microbiome_help.R"))

out_dir <- "/cwork/yx306/CorTree"
dmm_result_path <- file.path(out_dir, "REST_dmm_results.rds")

if (!file.exists(dmm_result_path)) {
  stop("Missing REST DMM result file: ", dmm_result_path)
}
dmm_result <- readRDS(dmm_result_path)

rest_paths <- list.files(
  out_dir,
  pattern = "^REST_task[0-9]+_.*_result_bundle\\.rds$",
  full.names = TRUE
)

if (length(rest_paths) == 0L) {
  stop("No REST result bundle files found in: ", out_dir)
}

rest_bundles <- lapply(rest_paths, readRDS) # all results, do not print

long_tables <- Map(
  function(bundle, path) {
    tab <- as.data.frame(bundle$ari_table, stringsAsFactors = FALSE)
    required_cols <- c("Method", "ARI")
    missing_cols <- setdiff(required_cols, names(tab))
    if (length(missing_cols) > 0L) {
      stop(
        "REST bundle ari_table is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        " in ",
        basename(path)
      )
    }
    out_tab <- data.frame(
      TaskID = rep.int(extract_task_id(path), nrow(tab)),
      InitMode = rep.int(bundle$init_mode, nrow(tab)),
      Method = as.character(tab$Method),
      MethodKey = ifelse(grepl("^Z_init", tab$Method), "Z_init", as.character(tab$Method)),
      ARI = as.numeric(tab$ARI),
      ResultFile = rep.int(basename(path), nrow(tab)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    dmm_ari <- adjusted_rand_index(
      bundle$data$chip_labels,
      dmm_result$cluster_assignments
    )

    rbind(
      out_tab,
      data.frame(
        TaskID = extract_task_id(path),
        InitMode = bundle$init_mode,
        Method = "DMM",
        MethodKey = "DMM",
        ARI = as.numeric(dmm_ari),
        ResultFile = basename(path),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  },
  rest_bundles,
  rest_paths
)

overall_long <- do.call(rbind.data.frame, c(long_tables, list(stringsAsFactors = FALSE)))
overall_long <- as.data.frame(overall_long, stringsAsFactors = FALSE, check.names = FALSE)
overall_long <- overall_long[order(overall_long[["TaskID"]], overall_long[["MethodKey"]]), , drop = FALSE]
row.names(overall_long) <- NULL

time_tables <- Map(
  function(bundle, path) {
    method_keys <- c("DMM", "IndTree", "CorTree")
    data.frame(
      TaskID = rep.int(extract_task_id(path), length(method_keys)),
      InitMode = rep.int(bundle$init_mode, length(method_keys)),
      MethodKey = method_keys,
      ElapsedSec = c(
        as.numeric(dmm_result$elapsed_time),
        vapply(
          c("IndTree", "CorTree"),
          function(method_key) extract_elapsed(bundle, method_key, bundle$init_mode),
          numeric(1)
        )
      ),
      ResultFile = rep.int(basename(path), length(method_keys)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  },
  rest_bundles,
  rest_paths
)

overall_time_long <- do.call(rbind.data.frame, c(time_tables, list(stringsAsFactors = FALSE)))
overall_time_long <- as.data.frame(overall_time_long, stringsAsFactors = FALSE, check.names = FALSE)
overall_time_long <- overall_time_long[
  order(overall_time_long[["TaskID"]], overall_time_long[["MethodKey"]]),
  ,
  drop = FALSE
]
row.names(overall_time_long) <- NULL

wide_input <- overall_long[, c("MethodKey", "InitMode", "ARI")]
if (nrow(wide_input) == 0L) {
  stop("No ARI rows found for REST bundles after filtering.")
}

time_wide_input <- overall_time_long[, c("MethodKey", "InitMode", "ElapsedSec")]
if (nrow(time_wide_input) == 0L) {
  stop("No elapsed-time rows found for REST bundles after filtering.")
}

duplicate_pairs <- duplicated(wide_input[c("MethodKey", "InitMode")])
if (any(duplicate_pairs)) {
  dup_keys <- unique(
    paste(
      wide_input$MethodKey[duplicate_pairs],
      wide_input$InitMode[duplicate_pairs],
      sep = " @ "
    )
  )
  stop(
    "Duplicate MethodKey/InitMode combinations found in REST results: ",
    paste(dup_keys, collapse = ", ")
  )
}

wide_matrix <- with(
  wide_input,
  tapply(ARI, list(MethodKey, InitMode), identity)
)
overall_wide_raw <- data.frame(
  MethodKey = row.names(wide_matrix),
  as.data.frame.matrix(wide_matrix),
  row.names = NULL,
  check.names = FALSE
)

init_mode_order <- unique(overall_long$InitMode[order(overall_long$TaskID)])
wide_cols <- intersect(init_mode_order, names(overall_wide_raw))
overall_wide_raw <- overall_wide_raw[, c("MethodKey", wide_cols), drop = FALSE]
row.names(overall_wide_raw) <- NULL

ari_matrix <- as.matrix(overall_wide_raw[, wide_cols, drop = FALSE])
storage.mode(ari_matrix) <- "numeric"
mean_ari <- rowMeans(ari_matrix, na.rm = TRUE)
nonmissing_counts <- rowSums(!is.na(ari_matrix))
constant_rows <- rowSums(
  abs(sweep(ari_matrix, 1L, mean_ari, FUN = "-")) > 1e-12,
  na.rm = TRUE
) == 0L & nonmissing_counts > 0L

common_ari <- rep(NA_real_, nrow(overall_wide_raw))
common_ari[constant_rows] <- mean_ari[constant_rows]
ari_matrix_display <- ari_matrix
ari_matrix_display[constant_rows, ] <- NA_real_

overall_wide <- data.frame(
  MethodKey = overall_wide_raw$MethodKey,
  CommonARI = common_ari,
  as.data.frame(ari_matrix_display, check.names = FALSE),
  MeanARI = mean_ari,
  row.names = NULL,
  check.names = FALSE
)

overall_wide <- overall_wide[order(-overall_wide$MeanARI), , drop = FALSE]
row.names(overall_wide) <- NULL

overall_wide_display <- overall_wide
overall_wide_display[] <- lapply(overall_wide_display, function(col) {
  if (is.numeric(col)) {
    replace(col, is.na(col), "")
  } else {
    replace(as.character(col), is.na(col), "")
  }
})

time_wide_matrix <- with(
  time_wide_input,
  tapply(ElapsedSec, list(MethodKey, InitMode), identity)
)
overall_time_wide_raw <- data.frame(
  MethodKey = row.names(time_wide_matrix),
  as.data.frame.matrix(time_wide_matrix),
  row.names = NULL,
  check.names = FALSE
)
overall_time_wide_raw <- overall_time_wide_raw[, c("MethodKey", wide_cols), drop = FALSE]
row.names(overall_time_wide_raw) <- NULL

elapsed_matrix <- as.matrix(overall_time_wide_raw[, wide_cols, drop = FALSE])
storage.mode(elapsed_matrix) <- "numeric"
mean_elapsed <- rowMeans(elapsed_matrix, na.rm = TRUE)
nonmissing_time_counts <- rowSums(!is.na(elapsed_matrix))
constant_time_rows <- rowSums(
  abs(sweep(elapsed_matrix, 1L, mean_elapsed, FUN = "-")) > 1e-12,
  na.rm = TRUE
) == 0L & nonmissing_time_counts > 0L

common_elapsed <- rep(NA_real_, nrow(overall_time_wide_raw))
common_elapsed[constant_time_rows] <- mean_elapsed[constant_time_rows]
elapsed_matrix_display <- elapsed_matrix
elapsed_matrix_display[constant_time_rows, ] <- NA_real_

overall_time_wide <- data.frame(
  MethodKey = overall_time_wide_raw$MethodKey,
  CommonElapsedSec = common_elapsed,
  as.data.frame(elapsed_matrix_display, check.names = FALSE),
  MeanElapsedSec = mean_elapsed,
  row.names = NULL,
  check.names = FALSE
)
overall_time_wide <- overall_time_wide[order(overall_time_wide$MethodKey), , drop = FALSE]
row.names(overall_time_wide) <- NULL

overall_time_wide_display <- overall_time_wide
overall_time_wide_display[] <- lapply(overall_time_wide_display, function(col) {
  if (is.numeric(col)) {
    replace(col, is.na(col), "")
  } else {
    replace(as.character(col), is.na(col), "")
  }
})

long_out <- file.path(out_dir, "REST_overall_ari_long.csv")
wide_out <- file.path(out_dir, "REST_overall_ari_table.csv")
time_long_out <- file.path(out_dir, "REST_overall_elapsed_long.csv")
time_wide_out <- file.path(out_dir, "REST_overall_elapsed_table.csv")
summary_out <- file.path(out_dir, "REST_overall_ari_summary.txt")

utils::write.csv(overall_long, long_out, row.names = FALSE)
utils::write.csv(overall_wide, wide_out, row.names = FALSE)
utils::write.csv(overall_time_long, time_long_out, row.names = FALSE)
utils::write.csv(overall_time_wide, time_wide_out, row.names = FALSE)
writeLines(
  c(
    "REST overall ARI table over different initial values",
    "",
    capture.output(print(overall_wide_display, row.names = FALSE)),
    "",
    "REST elapsed-time table (seconds; DMM, CorTree, and IndTree)",
    "",
    capture.output(print(overall_time_wide_display, row.names = FALSE)),
    "",
    "Detailed long-format table:",
    capture.output(print(overall_long, row.names = FALSE)),
    "",
    "Detailed elapsed-time long-format table (seconds; DMM, CorTree, and IndTree):",
    capture.output(print(overall_time_long, row.names = FALSE))
  ),
  con = summary_out
)

cat("REST analysis complete\n")
cat("Saved long-format table:", long_out, "\n")
cat("Saved overall ARI table:", wide_out, "\n")
cat("Saved elapsed-time long-format table:", time_long_out, "\n")
cat("Saved elapsed-time table:", time_wide_out, "\n")
cat("Saved summary text:", summary_out, "\n")
print(overall_wide_display)
print(overall_time_wide_display)

# ----------------------------------------------------------
# reorganize cluster-mean and ChIP-count plots across tasks
# ----------------------------------------------------------
library(ggplot2)
library(patchwork)

rest_plot_paths <- sub("_result_bundle\\.rds$", "_plot_artifacts.rds", rest_paths)
missing_plot_paths <- rest_plot_paths[!file.exists(rest_plot_paths)]
if (length(missing_plot_paths) > 0L) {
  stop(
    "Missing REST plot artifact files: ",
    paste(basename(missing_plot_paths), collapse = ", ")
  )
}

rest_task_info <- data.frame(
  TaskID = vapply(rest_paths, extract_task_id, integer(1)),
  InitMode = vapply(rest_bundles, function(bundle) bundle$init_mode, character(1)),
  PlotPath = rest_plot_paths,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rest_task_info <- rest_task_info[order(rest_task_info$TaskID), , drop = FALSE]
row.names(rest_task_info) <- NULL

resize_plot_text <- function(plot_obj) {
  plot_obj +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.key.size = grid::unit(0.65, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

rename_plot_title <- function(plot_obj, title_text) {
  plot_obj + labs(title = title_text)
}

fixed_method_order <- c("chip_labels", "kmeans", "pam", "centipede")
fixed_method_titles <- c(
  chip_labels = "ChIP_label",
  kmeans = "Kmeans",
  pam = "PAM",
  centipede = "CENTIPEDE"
)

task_specific_order <- c("init", "indtree", "cortree")
task_specific_titles <- c(
  init = "Z_init",
  indtree = "IndTree",
  cortree = "CorTree"
)

base_artifacts <- readRDS(rest_task_info$PlotPath[[1L]])
cluster_mean_panels <- lapply(
  fixed_method_order,
  function(method_key) {
    rename_plot_title(
      resize_plot_text(base_artifacts$cluster_mean[[method_key]]),
      fixed_method_titles[[method_key]]
    )
  }
)
histogram_panels <- lapply(
  fixed_method_order,
  function(method_key) {
    rename_plot_title(
      resize_plot_text(base_artifacts$histogram[[method_key]]),
      fixed_method_titles[[method_key]]
    )
  }
)

combined_plot_list <- c(cluster_mean_panels, histogram_panels)
rm(base_artifacts)
gc()

for (idx in seq_len(nrow(rest_task_info))) {
  task_label <- paste0("Task ", rest_task_info$TaskID[[idx]], " (", rest_task_info$InitMode[[idx]], ")")
  task_artifacts <- readRDS(rest_task_info$PlotPath[[idx]])

  task_cluster_row <- lapply(
    task_specific_order,
    function(method_key) {
      rename_plot_title(
        resize_plot_text(task_artifacts$cluster_mean[[method_key]]),
        paste(task_label, task_specific_titles[[method_key]], sep = "\n")
      )
    }
  )
  task_hist_row <- lapply(
    task_specific_order,
    function(method_key) {
      rename_plot_title(
        resize_plot_text(task_artifacts$histogram[[method_key]]),
        paste(task_label, task_specific_titles[[method_key]], sep = "\n")
      )
    }
  )

  combined_plot_list <- c(
    combined_plot_list,
    task_cluster_row,
    list(plot_spacer()),
    task_hist_row,
    list(plot_spacer())
  )

  rm(task_artifacts)
  gc()
}

rest_reorganized_plot <- wrap_plots(combined_plot_list, ncol = 4)
rest_reorganized_plot_out <- file.path(out_dir, "REST_reorganized_method_comparison.pdf")
rest_reorganized_plot_out_png <- file.path(out_dir, "REST_reorganized_method_comparison.png")

ggplot2::ggsave(
  filename = rest_reorganized_plot_out,
  plot = rest_reorganized_plot,
  width = 18,
  height = 28
)
ggplot2::ggsave(
  filename = rest_reorganized_plot_out_png,
  plot = rest_reorganized_plot,
  width = 18,
  height = 28,
  dpi = 300
)

cat("Saved reorganized method-comparison plot:", rest_reorganized_plot_out, "\n")
cat("Saved reorganized method-comparison plot:", rest_reorganized_plot_out_png, "\n")

# ----------------------------------------------------------
# redraw REST task 02 (pwm_quantile) for publication layout
# ----------------------------------------------------------
pub_plot_text <- function(plot_obj) {
  plot_obj +
    theme(
      plot.title = element_text(size = 15, hjust = 0.5),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.key.size = grid::unit(0.55, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

task02_idx <- which(
  rest_task_info$TaskID == 2L &
    rest_task_info$InitMode == "pwm_quantile"
)
if (length(task02_idx) != 1L) {
  stop("Could not uniquely identify REST task 02 pwm_quantile plot artifacts.")
}

task02_artifacts <- readRDS(rest_task_info$PlotPath[[task02_idx]])
task02_bundle <- rest_bundles[[match(rest_task_info$TaskID[[task02_idx]], vapply(rest_paths, extract_task_id, integer(1)))]]

dmm_cluster_plot <- rename_plot_title(
  pub_plot_text(
    cluster_mean_gg(
      task02_bundle$data$X,
      dmm_result$cluster_assignments,
      y_max = NULL
    )
  ),
  "DMM"
)
dmm_hist_plot <- rename_plot_title(
  pub_plot_text(
    hist_by_group_gg(
      task02_bundle$data$chip_data,
      dmm_result$cluster_assignments,
      title = NULL,
      y_max = 800
    ) + labs(fill = NULL)
  ),
  "DMM"
)

task02_cluster_row_fixed <- lapply(
  fixed_method_order,
  function(method_key) {
    rename_plot_title(
      pub_plot_text(task02_artifacts$cluster_mean[[method_key]]),
      fixed_method_titles[[method_key]]
    )
  }
)
task02_hist_row_fixed <- lapply(
  fixed_method_order,
  function(method_key) {
    rename_plot_title(
      pub_plot_text(task02_artifacts$histogram[[method_key]] + labs(fill = NULL)),
      fixed_method_titles[[method_key]]
    )
  }
)

task02_specific_order <- c("init", "indtree", "cortree")
task02_specific_titles <- c(
  init = "Z_init",
  indtree = "IndTree",
  cortree = "CorTree"
)

task02_cluster_row_model <- c(
  list(dmm_cluster_plot),
  lapply(
    task02_specific_order,
    function(method_key) {
      rename_plot_title(
        pub_plot_text(task02_artifacts$cluster_mean[[method_key]]),
        task02_specific_titles[[method_key]]
      )
    }
  )
)
task02_hist_row_model <- c(
  list(dmm_hist_plot),
  lapply(
    task02_specific_order,
    function(method_key) {
      rename_plot_title(
        pub_plot_text(task02_artifacts$histogram[[method_key]] + labs(fill = NULL)),
        task02_specific_titles[[method_key]]
      )
    }
  )
)

task02_pub_plot <- wrap_plots(
  c(
    task02_cluster_row_fixed,
    task02_hist_row_fixed,
    task02_cluster_row_model,
    task02_hist_row_model
  ),
  ncol = 4
)

task02_pub_plot_out <- file.path(out_dir, "REST_task02_pwm_quantile_reorganized.pdf")

ggplot2::ggsave(
  filename = task02_pub_plot_out,
  plot = task02_pub_plot,
  width = 16,
  height = 16
)
task02_pub_plot_out_png <- file.path(out_dir, "REST_task02_pwm_quantile_reorganized.png")
ggplot2::ggsave(
  filename = task02_pub_plot_out_png,
  plot = task02_pub_plot,
  width = 16,
  height = 16,
  dpi = 300
)

cat("Saved REST task 02 reorganized plot:", task02_pub_plot_out, "\n")
cat("Saved REST task 02 reorganized plot:", task02_pub_plot_out_png, "\n")


# ----------------------------------------------------------
# check convergence for different initial values for CorTree
# ----------------------------------------------------------
## Full CorTree mixing-weight traceplots with ggplot
## One row of subplots, one panel per initialization mode

library(tidyr)
library(dplyr)
ordered_idx <- order(vapply(rest_paths, extract_task_id, integer(1)))
rest_bundles_ordered <- rest_bundles[ordered_idx] # all data, do not print

trace_plot_list <- lapply(rest_bundles_ordered, function(bundle) {
  pi_chain <- as.matrix(bundle$fit$cortree$mcmc$pi)
  storage.mode(pi_chain) <- "numeric"

  if (is.null(bundle$data$chip_labels) || is.null(bundle$init$selected) || is.null(bundle$cluster_assignments$cortree)) {
    stop("REST bundle is missing chip_labels, selected initialization labels, or CorTree cluster assignments.")
  }
  init_ari <- adjusted_rand_index(bundle$data$chip_labels, bundle$init$selected)
  cortree_ari <- adjusted_rand_index(bundle$data$chip_labels, bundle$cluster_assignments$cortree)
  subplot_title <- sprintf(
    "Init: %s\nARI(init, ChIP): %.3f\nARI(CorTree, ChIP): %.3f",
    bundle$init_mode,
    init_ari,
    cortree_ari
  )

  trace_df <- as.data.frame(t(pi_chain), check.names = FALSE)
  names(trace_df) <- paste0("Cluster ", seq_len(ncol(trace_df)))
  trace_df$Iteration <- seq_len(nrow(trace_df))

  trace_long <- tidyr::pivot_longer(
    trace_df,
    cols = starts_with("Cluster "),
    names_to = "Cluster",
    values_to = "Value"
  )

  ggplot(trace_long, aes(x = Iteration, y = Value, color = Cluster)) +
    geom_line(linewidth = 0.45) +
    geom_vline(xintercept = 100, linetype = "dashed", color = "grey40") +
    scale_color_manual(
      values = c(
        "Cluster 1" = "#F8766D",
        "Cluster 2" = "#A3A500",
        "Cluster 3" = "#00BF7D",
        "Cluster 4" = "#00B0F6",
        "Cluster 5" = "#E76BF3"
      )
    ) +
    labs(
      title = subplot_title,
      x = "Iteration",
      y = "Value"
    ) +
    scale_y_continuous(limits = c(0, 1))+
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
})

combined_trace_plot <- wrap_plots(trace_plot_list, nrow = 1, guides = "collect") +
  theme(legend.position = "right")

trace_plot_out <- file.path(out_dir, "REST_cortree_pi_trace_full_ggplot.pdf")
ggplot2::ggsave(
  filename = trace_plot_out,
  plot = combined_trace_plot,
  width = 14,
  height = 4.2
)

trace_plot_out_png <- file.path(out_dir, "REST_cortree_pi_trace_full_ggplot.png")
ggplot2::ggsave(
  filename = trace_plot_out_png,
  plot = combined_trace_plot,
  width = 14,
  height = 4.2,
  dpi = 300
)

cat("Saved full CorTree ggplot traceplots:", trace_plot_out, "\n")
cat("Saved full CorTree ggplot traceplots:", trace_plot_out_png, "\n")
