plot_pos_neg_profiles <- function(X, chip_labels, y_max = NULL,n_thresh = 10, title = "") {
  chip_labels = round(chip_labels)
  # Convert labels to character (so we can index table easily)
  lab_char <- as.character(chip_labels)
  counts   <- table(lab_char)
  
  # Keep only those labels whose count exceeds n_thresh
  valid_labels <- names(counts)[counts > n_thresh]
  if (length(valid_labels) == 0) {
    stop("No labels exceed the threshold of ", n_thresh)
  }
  
  # Compute one mean‐profile per valid label
  profiles <- sapply(valid_labels, function(lbl) {
    colMeans(X[lab_char == lbl, , drop = FALSE], na.rm = TRUE)
  })
  # ‘profiles’ is a matrix: ncol(X) rows × length(valid_labels) columns
  
  x <- seq_len(nrow(profiles))
  y_min <- 0
  
  if(is.null(y_max)){
    y_max <- max(profiles, na.rm = TRUE)
  }
  
  cols <- rainbow(length(valid_labels))
  # Plot the first profile
  plot(
    x, profiles[, 1],
    type  = "l",
    col   = adjustcolor(cols[1], alpha.f = 0.7),  # add 70% opacity here
    lwd   = 2,
    ylim  = c(y_min, y_max),
    xlab  = "Position",
    ylab  = "Signal",
    main  = title
  )
  
  # Add the remaining profiles
  if (length(valid_labels) > 1) {
    for (i in 2:length(valid_labels)) {
      lines(
        x, profiles[, i],
        col = adjustcolor(cols[i], alpha.f = 0.7),  # add 70% opacity here
        lwd = 2
      )
    }
  }
  
  # Legend entries are “label:count”
  legend_text <- paste0(valid_labels, ":", counts[valid_labels])
  legend(
    "topright",
    legend = legend_text,
    col    = cols,
    lty    = 1,
    lwd    = 2,
    bty    = "n"
  )
}

hist_by_group_gg <- function(x, z,
                             bins     = 30,
                             palette  = "Set1",
                             alpha    = 0.5,
                             n_thresh = 10,
                             y_max    = NULL,
                             title    = NULL,
                             xlab     = NULL,
                             ylab     = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the ggplot2 package to use this function.")
  }
  
  # Build initial data.frame
  df <- data.frame(
    value = log(x + 1),
    group = factor(z)
  )
  
  # Filter groups by threshold
  counts    <- table(df$group)
  valid_grp <- names(counts)[counts >= n_thresh]
  if (length(valid_grp) == 0) {
    stop("No groups have at least ", n_thresh, " observations.")
  }
  df <- subset(df, group %in% valid_grp)
  df$group <- droplevels(df$group)
  
  # Recompute legend labels
  grp_levels   <- levels(df$group)
  grp_counts   <- as.vector(table(df$group)[grp_levels])
  legend_labels <- paste0(grp_levels, ":", grp_counts)
  
  # Base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) +
    ggplot2::geom_histogram(position = "identity",
                            bins     = bins,
                            alpha    = alpha) +
    ggplot2::scale_fill_brewer(palette = palette,
                               labels  = legend_labels) +
    ggplot2::labs(
      title = if (is.null(title)) "Histogram by Group" else title,
      x     = if (is.null(xlab))   "log(x+1)" else xlab,
      y     = if (is.null(ylab))   "Count"    else ylab,
      fill  = "Group"
    ) +
    ggplot2::theme_minimal()+
    ggplot2::theme(
      legend.position = "top",
      legend.title    = ggplot2::element_blank()
    )
  
  # Apply y-axis limit if requested
  if (!is.null(y_max)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, y_max))
  }
  
  return(p)
}


cluster_mean_gg <- function(X, chip_labels,
                            y_max    = NULL,
                            n_thresh = 10,
                            title    = "") {
  # Dependencies
  if (!requireNamespace("ggplot2", quietly=TRUE) ||
      !requireNamespace("tidyr",   quietly=TRUE)) {
    stop("Please install the ggplot2 and tidyr packages to use this function.")
  }
  
  # 1) Round & tabulate
  lab_char <- as.character(round(chip_labels))
  counts   <- table(lab_char)
  
  # 2) Filter by threshold
  valid <- names(counts)[counts > n_thresh]
  if (length(valid) == 0) {
    stop("No labels exceed the threshold of ", n_thresh)
  }
  
  # 3) Compute mean‐profiles (columns = positions)
  profiles <- sapply(valid, function(lbl) {
    colMeans(X[lab_char == lbl, , drop=FALSE], na.rm=TRUE)
  })
  # profiles: ncol(X) rows × length(valid) columns
  
  # 4) Build long data frame
  df <- as.data.frame(profiles)
  df$Position <- seq_len(nrow(df))
  df_long <- tidyr::pivot_longer(
    df,
    cols      = -Position,
    names_to  = "label",
    values_to = "signal"
  )
  
  # 5) Factor labels to include counts
  df_long$label <- factor(
    df_long$label,
    levels = valid,
    labels = paste0(valid, ":", counts[valid])
  )
  
  # 6) Auto‐compute y_max if needed
  if (is.null(y_max)) {
    y_max <- max(df_long$signal, na.rm=TRUE)
  }
  
  # 7) Plot
  p <- ggplot2::ggplot(df_long,
                       ggplot2::aes(x = Position,
                                    y = signal,
                                    color = label)) +
    ggplot2::geom_line(linewidth = 1, alpha = 0.7) +
    ggplot2::scale_color_brewer(palette = "Set1", name = NULL) +
    ggplot2::labs(
      title = title,
      x     = "Position",
      y     = "Signal"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_max)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      legend.title    = ggplot2::element_blank()
    )
  
  return(p)
}
