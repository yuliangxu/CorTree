plot_pos_neg_profiles <- function(pos, neg, title = "") {
  # pos, neg: numeric vectors of equal length
  if (length(pos) != length(neg)) {
    stop("‘pos’ and ‘neg’ must be the same length.")
  }
  
  x <- seq_along(pos)
  y_min <- min(c(pos, neg))
  y_max <- max(c(pos, neg))
  
  # First, plot the positive profile
  plot(
    x, pos,
    type = "l",
    col  = "blue",
    lwd  = 2,
    ylim = c(y_min, y_max),
    xlab = "Position",
    ylab = "Signal",
    main = title
  )
  
  # Then add the negative profile
  lines(
    x, neg,
    col = "red",
    lwd = 2
  )
  
  # Add a legend in the top-right corner
  legend(
    "topright",
    legend = c("Positive", "Negative"),
    col    = c("blue", "red"),
    lty    = 1,
    lwd    = 2,
    bty    = "n"
  )
}