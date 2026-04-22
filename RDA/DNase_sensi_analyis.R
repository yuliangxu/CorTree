rm(list = ls())

project_dir <- normalizePath(".")
source(file.path(project_dir, "R", "DCC_help.R"))
source(file.path(project_dir, "R", "microbiome_help.R"))

library(ggplot2)
library(patchwork)

default_out_dir <- "/cwork/yx306/CorTree"
summary_out_dir <- default_out_dir
dir.create(summary_out_dir, recursive = TRUE, showWarnings = FALSE)

sensitivity_label_map <- c(
  cutoff_layer_4 = "cutoff_layer=4",
  c_sigma2_vec_1 = "c_sigma2_vec=1",
  sigma_mu2_1 = "sigma_mu2=1"
)

rest_records <- collect_sensitivity_dataset_results(
  "REST",
  default_out_dir,
  default_out_dir,
  sensitivity_label_map
)
nrf1_records <- collect_sensitivity_dataset_results(
  "NRF1",
  default_out_dir,
  default_out_dir,
  sensitivity_label_map
)
rest_records <- filter_sensitivity_records(rest_records)
nrf1_records <- filter_sensitivity_records(nrf1_records)
combined_records <- rbind(rest_records, nrf1_records)

rest_task_summary <- build_sensitivity_task_summary(rest_records)
nrf1_task_summary <- build_sensitivity_task_summary(nrf1_records)
combined_task_summary <- rbind(rest_task_summary, nrf1_task_summary)

rest_method_wide <- build_sensitivity_method_wide(rest_records)
nrf1_method_wide <- build_sensitivity_method_wide(nrf1_records)

combined_compare <- merge(
  rest_task_summary[, c("SensitivityMode", "SensitivityLabel", "CorTree_ARI", "ARI_Improvement"), drop = FALSE],
  nrf1_task_summary[, c("SensitivityMode", "SensitivityLabel", "CorTree_ARI", "ARI_Improvement"), drop = FALSE],
  by = c("SensitivityMode", "SensitivityLabel"),
  suffixes = c("_REST", "_NRF1"),
  all = TRUE,
  sort = FALSE
)

combined_long_out <- file.path(summary_out_dir, "DNase_sensitivity_ari_long.csv")
combined_task_out <- file.path(summary_out_dir, "DNase_sensitivity_task_summary.csv")
rest_wide_out <- file.path(summary_out_dir, "REST_sensitivity_ari_table.csv")
nrf1_wide_out <- file.path(summary_out_dir, "NRF1_sensitivity_ari_table.csv")
combined_compare_out <- file.path(summary_out_dir, "DNase_sensitivity_dataset_comparison.csv")
summary_txt_out <- file.path(summary_out_dir, "DNase_sensitivity_summary.txt")
rest_plot_out <- file.path(summary_out_dir, "REST_sensitivity_combined_plots.pdf")
rest_plot_out_png <- file.path(summary_out_dir, "REST_sensitivity_combined_plots.png")
nrf1_plot_out <- file.path(summary_out_dir, "NRF1_sensitivity_combined_plots.pdf")
nrf1_plot_out_png <- file.path(summary_out_dir, "NRF1_sensitivity_combined_plots.png")

utils::write.csv(combined_records, combined_long_out, row.names = FALSE)
utils::write.csv(combined_task_summary, combined_task_out, row.names = FALSE)
utils::write.csv(rest_method_wide, rest_wide_out, row.names = FALSE)
utils::write.csv(nrf1_method_wide, nrf1_wide_out, row.names = FALSE)
utils::write.csv(combined_compare, combined_compare_out, row.names = FALSE)

rest_combined_plot <- build_combined_sensitivity_dataset_plot("REST", rest_records, default_out_dir)
nrf1_combined_plot <- build_combined_sensitivity_dataset_plot("NRF1", nrf1_records, default_out_dir)

ggplot2::ggsave(filename = rest_plot_out, plot = rest_combined_plot, width = 18, height = 8.5)
ggplot2::ggsave(filename = rest_plot_out_png, plot = rest_combined_plot, width = 18, height = 8.5, dpi = 300)
ggplot2::ggsave(filename = nrf1_plot_out, plot = nrf1_combined_plot, width = 18, height = 8.5)
ggplot2::ggsave(filename = nrf1_plot_out_png, plot = nrf1_combined_plot, width = 18, height = 8.5, dpi = 300)

writeLines(
  c(
    "DNase sensitivity analysis ARI summary",
    "",
    format_sensitivity_dataset_summary(rest_task_summary),
    format_sensitivity_dataset_summary(nrf1_task_summary),
    "REST task summary:",
    capture.output(print(rest_task_summary, row.names = FALSE)),
    "",
    "NRF1 task summary:",
    capture.output(print(nrf1_task_summary, row.names = FALSE)),
    "",
    "REST method-by-sensitivity ARI table:",
    capture.output(print(rest_method_wide, row.names = FALSE)),
    "",
    "NRF1 method-by-sensitivity ARI table:",
    capture.output(print(nrf1_method_wide, row.names = FALSE)),
    "",
    "Cross-dataset CorTree comparison:",
    capture.output(print(combined_compare, row.names = FALSE))
  ),
  con = summary_txt_out
)

cat("DNase sensitivity analysis summary complete\n")
cat("Saved combined long-format table:", combined_long_out, "\n")
cat("Saved task summary table:", combined_task_out, "\n")
cat("Saved REST sensitivity ARI table:", rest_wide_out, "\n")
cat("Saved NRF1 sensitivity ARI table:", nrf1_wide_out, "\n")
cat("Saved cross-dataset comparison:", combined_compare_out, "\n")
cat("Saved REST combined plot:", rest_plot_out, "\n")
cat("Saved REST combined plot:", rest_plot_out_png, "\n")
cat("Saved NRF1 combined plot:", nrf1_plot_out, "\n")
cat("Saved NRF1 combined plot:", nrf1_plot_out_png, "\n")
cat("Saved text summary:", summary_txt_out, "\n")
print(rest_task_summary)
print(nrf1_task_summary)
