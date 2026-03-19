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
