construct_tree <- function(X, tree_depth) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  if (length(tree_depth) != 1L || !is.finite(tree_depth) || tree_depth < 1) {
    stop("`tree_depth` must be a positive scalar.")
  }

  m <- as.integer(tree_depth)
  n <- nrow(X)
  k <- ncol(X)
  total_nodes <- 2L^(m + 1L) - 1L
  total_parents <- 2L^m - 1L

  count <- matrix(0, nrow = n, ncol = total_nodes)
  kappa <- matrix(0, nrow = n, ncol = total_parents)
  interval_left <- integer(total_nodes)
  interval_right <- integer(total_nodes)

  count[, 1L] <- rowSums(X)
  interval_left[1L] <- 0L
  interval_right[1L] <- as.integer(k)

  for (node in seq_len(total_parents)) {
    left_child <- 2L * node
    right_child <- left_child + 1L

    interval_left[left_child] <- interval_left[node]
    interval_right[left_child] <- (interval_left[node] + interval_right[node]) %/% 2L
    interval_left[right_child] <- interval_right[left_child]
    interval_right[right_child] <- interval_right[node]

    if (interval_left[node] < interval_right[left_child]) {
      left_idx <- seq.int(interval_left[node] + 1L, interval_right[left_child])
      count[, left_child] <- rowSums(X[, left_idx, drop = FALSE])
    }

    if (interval_left[right_child] < interval_right[node]) {
      right_idx <- seq.int(interval_left[right_child] + 1L, interval_right[node])
      count[, right_child] <- rowSums(X[, right_idx, drop = FALSE])
    }

    kappa[, node] <- count[, left_child] - count[, node] / 2
  }

  parent_count <- count[, seq_len(total_parents), drop = FALSE]
  left_count <- count[, 2L * seq_len(total_parents), drop = FALSE]
  empirical_phi <- qlogis((left_count + 0.5) / (parent_count + 1.0))

  list(
    count = count,
    kappa = kappa,
    interval_left = interval_left,
    interval_right = interval_right,
    empirical_phi = empirical_phi
  )
}
