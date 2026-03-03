rm(list = ls())

sim_dir <- "./sim_result"
files <- list.files(sim_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0L) {
  stop("No .rds files found in ./sim_result")
}

# Expected column order in saved matrices:
# 1: K-means, 2: PAM, 3: IndTree, 4: CorTree
method_names <- c("K-means", "PAM", "IndTree", "CorTree")

make_summary_table <- function(ari_mat, time_mat, method_names, file_name = NA_character_) {
  out <- data.frame(
    File = file_name,
    Method = method_names,
    ARI_mean = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, mean, na.rm = TRUE),
    ARI_sd = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, sd, na.rm = TRUE),
    ARI_median = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, median, na.rm = TRUE),
    ARI_n_le_0 = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, function(x) sum(x <= 0, na.rm = TRUE)),
    ARI_n_eq_1 = apply(ari_mat[, seq_along(method_names), drop = FALSE], 2, function(x) sum(x == 1, na.rm = TRUE)),
    Avg_time_sec = apply(time_mat[, seq_along(method_names), drop = FALSE], 2, mean, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out$ARI_mean <- round(out$ARI_mean, 4)
  out$ARI_sd <- round(out$ARI_sd, 4)
  out$ARI_median <- round(out$ARI_median, 4)
  out$Avg_time_sec <- round(out$Avg_time_sec, 3)
  out
}

all_ari <- vector("list", length(method_names))
all_time <- vector("list", length(method_names))
per_file_tables <- list()

for (i in seq_along(method_names)) {
  all_ari[[i]] <- numeric(0)
  all_time[[i]] <- numeric(0)
}

for (f in files) {
  obj <- readRDS(f)
  if (!is.list(obj) || !all(c("ari", "time") %in% names(obj))) {
    warning("Skipping file with unexpected structure: ", basename(f))
    next
  }

  ari_mat <- as.matrix(obj$ari)
  time_mat <- as.matrix(obj$time)

  if (ncol(ari_mat) < length(method_names) || ncol(time_mat) < length(method_names)) {
    warning("Skipping file with fewer than 4 method columns: ", basename(f))
    next
  }

  # Subtable for this file
  one_tbl <- make_summary_table(
    ari_mat = ari_mat,
    time_mat = time_mat,
    method_names = method_names,
    file_name = basename(f)
  )
  per_file_tables[[basename(f)]] <- one_tbl

  cat("\nSubtable for", basename(f), "\n")
  print(one_tbl[, -1, drop = FALSE], row.names = FALSE)

  out_one_csv <- file.path(
    sim_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_summary_table.csv")
  )
  write.csv(one_tbl, out_one_csv, row.names = FALSE)

  for (j in seq_along(method_names)) {
    all_ari[[j]] <- c(all_ari[[j]], ari_mat[, j])
    all_time[[j]] <- c(all_time[[j]], time_mat[, j])
  }
}

summary_tbl <- make_summary_table(
  ari_mat = do.call(cbind, all_ari),
  time_mat = do.call(cbind, all_time),
  method_names = method_names,
  file_name = "ALL_FILES"
)

cat("\nFinal summary across all files in ./sim_result\n")
print(summary_tbl[, -1, drop = FALSE], row.names = FALSE)

out_csv <- file.path(sim_dir, "Analyze1_sim1_final_table.csv")
write.csv(summary_tbl, out_csv, row.names = FALSE)
cat("\nSaved:", out_csv, "\n")
