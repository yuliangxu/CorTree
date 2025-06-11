# Required Libraries
library(igraph)
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(gridExtra)
# with 95% quantiles
# Load required libraries
library(dplyr)
library(tidyr)

# Function to check if matrix is symmetric positive definite
isSymmetricPositiveDefinite <- function(A) {
  # Check if the matrix is symmetric
  sym_bool = isSymmetric(A)
  eigenvalues <- eigen(A, only.values = TRUE)$values
  pd_bool = all(eigenvalues > 0)
  
  print(paste("Symmetric: ", sym_bool,"; Positive definite: ", pd_bool))
  
}
# Function to build dyadic tree
plot_dyadic_tree <- function(X, vertex.size = 15, vertex.label.cex = 1) {
  # Number of nodes in the vector
  n <- length(X)
  
  # Create an empty graph
  g <- make_empty_graph(n = 0, directed = TRUE)
  
  # Add nodes (vertices) to the graph
  g <- add_vertices(g, n, label = as.character(X))
  
  # Add edges based on the parent-child relationships
  for (i in 1:n) {
    left_child <- 2 * (i-1) + 2
    right_child <- 2 * (i-1) + 3

    # Add an edge from the current node to its left child if exists
    if (left_child <= n) {
      g <- add_edges(g, c(i, left_child))
    }
    
    # Add an edge from the current node to its right child if exists
    if (right_child <= n) {
      g <- add_edges(g, c(i, right_child))
    }
  }
  
  # Define a tree layout with the 'layout.reingold.tilford' method for improved spacing
  layout <- layout_with_sugiyama(g)$layout
  # layout <- layout_as_tree(g)
  
  # Plot the graph with no overlapping of nodes
  plot(g, vertex.label = V(g)$label, layout = layout, 
       edge.arrow.size = 0.2, vertex.size = vertex.size, 
       vertex.label.cex = vertex.label.cex,
       asp = 0.6,
       vertex.color = "white" )
}



# generate toy example for count data 
# with n sample, k clutser, and p total nodes
gen_tree_fast = function(n,k,p,cutoff_layer) {
  total_nodes = p
  total_parents = floor((total_nodes - 1) / 2)
  total_correlated = 2^(cutoff_layer+1) - 1
  
  true_params$pi = rep(1/k,k) 
  true_params$pi = true_params$pi + 0.2*runif(k)
  true_params$pi = true_params$pi / sum(true_params$pi)
  total_correlated = 2^(cutoff_layer+1) - 1
  mu_1 = c(rep(-2,5),rep(2,5),rep(0,p-10))
  mu_2 = c(rep(0,5), rep(-2,5),rep(2,5),rep(0,p-15))
  true_params$mu = cbind(mu_1, mu_2)
  true_params$Sigma = array(0, dim = c(total_correlated,total_correlated,k))

  for(i in 1:k){
    row_nonzero = sample(1:total_correlated, floor(total_correlated*0.3))
    col_nonzero = row_nonzero
    Sigma_k = diag(total_correlated)
    Sigma_k[row_nonzero, col_nonzero] = rnorm(length(row_nonzero))
    true_params$Sigma[,,i] = t(Sigma_k) %*% Sigma_k
  }

  

  # generate one tree
  gen_one_tree = function(all_count,total_nodes,V){
    X = rep(0, total_nodes)
    total_parents = floor((total_nodes - 1) / 2)
    X[1] = all_count
    loglike = 0
    for (node in 1:total_parents) {
      # parent <- floor((i - 1) / 2) + 1
      left_child <- 2 * (node-1) + 2
      right_child <- 2 * (node-1) + 3

      X[left_child] = rbinom(1, X[node], V[node])
      X[right_child] = X[node] - X[left_child]
      
      loglike = loglike + 
        X[left_child]*log(V[node]) + X[right_child]*log(1-V[node])
    }
    
    return(list(sample = X, 
                loglike = loglike))
  }
  
  # generate tree structure
  count = matrix(0, nrow = n, ncol = p-total_parents)
  phi = matrix(0, nrow = n, ncol = p)
  tree_structure = matrix(0, nrow = n, ncol = p)
  loglike = 0
  
  Z = sample(1:k, n, prob = true_params$pi, replace=T)
  for(i_clus in 1:k){
    idx_k = which(Z==i_clus)
    n_k = sum(Z==i_clus)
    phi[idx_k,1:total_correlated] = MASS::mvrnorm(n_k, 
                            mu = true_params$mu[1:total_correlated,i_clus], 
                            Sigma = true_params$Sigma[1:total_correlated,1:total_correlated,i_clus])
    phi[idx_k,(total_correlated+1):total_nodes] = 
     t( t(matrix(rnorm(n_k*(total_nodes-total_correlated), sd = 0.1), 
             nrow = n_k, ncol = total_nodes-total_correlated)) + 
      true_params$mu[(total_correlated+1):total_nodes,i_clus] )
  }
  
  

  for(i in 1:n){
    Prob_i = 1/(1+exp(-phi[i,]))
    total_count_i = sample(1e3:1e4,1)
    tree_i = gen_one_tree(all_count = total_count_i,
                             total_nodes = total_nodes, V = Prob_i)
    tree_structure[i,] = tree_i$sample
    loglike = loglike + tree_i$loglike
    count[i,] = tree_structure[i,(total_parents+1):total_nodes]
  }
  
  true_params$phi = phi
  true_params$loglike = loglike
  true_params$Z = Z
  
  
  
  return(list(count = count, 
              tree_structure = tree_structure,
              true_params = true_params))
}


# generate toy example for count data 
# with n sample, k clutser, and p total nodes
gen_tree = function(n,k,p,cutoff_layer, true_params = NULL) {
  total_nodes = p
  total_parents = floor((total_nodes - 1) / 2)
  total_correlated = 2^(cutoff_layer+1) - 1
  
  if(is.null(true_params)) {
    true_params = list()
    true_params$pi = rep(1/k,k) 
    true_params$pi = true_params$pi + 0.1*runif(k)
    true_params$pi = true_params$pi / sum(true_params$pi)
    true_params$mu = matrix(runif(k*p), nrow = p, ncol = k)
    true_params$Sigma = array(0, dim = c(p,p,k))
    
    for(i in 1:k){
      row_nonzero = sample(1:total_correlated, floor(total_correlated*0.3))
      col_nonzero = row_nonzero
      Sigma_k = diag(p)
      Sigma_k[row_nonzero, col_nonzero] = rnorm(length(row_nonzero))
      true_params$Sigma[,,i] = t(Sigma_k) %*% Sigma_k
    }
  }
  
  
  
  # generate one tree
  gen_one_tree = function(all_count,total_nodes,V){
    X = rep(0, total_nodes)
    total_parents = floor((total_nodes - 1) / 2)
    X[1] = all_count
    loglike = 0
    for (node in 1:total_parents) {
      # parent <- floor((i - 1) / 2) + 1
      left_child <- 2 * (node-1) + 2
      right_child <- 2 * (node-1) + 3
      
      X[left_child] = rbinom(1, X[node], V[node])
      X[right_child] = X[node] - X[left_child]
      
      loglike = loglike + 
        X[left_child]*log(V[node]) + X[right_child]*log(1-V[node])
    }
    
    return(list(sample = X, 
                loglike = loglike))
  }
  
  # generate tree structure
  count = matrix(0, nrow = n, ncol = p-total_parents)
  phi = matrix(0, nrow = n, ncol = p)
  tree_structure = matrix(0, nrow = n, ncol = p)
  loglike = 0
  
  Z = sample(1:k, n, prob = true_params$pi, replace=T)
  
  
  for(i in 1:n){
    phi[i,] = MASS::mvrnorm(1,  
                            mu = true_params$mu[,Z[i]], 
                            Sigma = true_params$Sigma[,,Z[i]])
    Prob_i = 1/(1+exp(-phi[i,]))
    total_count_i = sample(1e3:1e4,1)
    tree_i = gen_one_tree(all_count = total_count_i,
                          total_nodes = total_nodes, V = Prob_i)
    tree_structure[i,] = tree_i$sample
    loglike = loglike + tree_i$loglike
    count[i,] = tree_structure[i,(total_parents+1):total_nodes]
  }
  
  true_params$phi = phi
  true_params$loglike = loglike
  true_params$Z = Z
  
  
  
  return(list(count = count, 
              tree_structure = tree_structure,
              true_params = true_params))
}


plot_two_mat = function(A,B){
  # Reshape matrices into data frames for ggplot
  A_melted <- melt(A)
  B_melted <- melt(B)
  
  # Add labels to distinguish the matrices
  A_melted$Matrix <- deparse(substitute(A))
  B_melted$Matrix <- deparse(substitute(B))
  
  # Combine both data frames
  data_combined <- rbind(A_melted, B_melted)
  
  # Get the global range for the color scale
  global_range <- range(data_combined$value)
  
  # Create the plot
  p <- ggplot(data_combined, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = global_range, name = "Value") +
    facet_wrap(~ Matrix) +
    scale_y_reverse() +  # Reverse the y-axis to match the matrix row order
    labs(x = "Column", y = "Row", 
         title = "Comparison of Symmetric Positive Definite Matrices") +
    theme_minimal()
  
  # Print the plot
  print(p)
}

plot_mat = function(A,title = NULL){
  # Reshape matrices into data frames for ggplot
  A_melted <- melt(A)
  
  # Add labels to distinguish the matrices
  A_melted$Matrix <- deparse(substitute(A))
  
  # Combine both data frames
  data_combined <- rbind(A_melted)
  
  # Get the global range for the color scale
  global_range <- range(data_combined$value)
  if(is.null(title)){
    title = deparse(substitute(A))
  }
  
  # Create the plot
  p <- ggplot(data_combined, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = global_range, name = "Value") +
    facet_wrap(~ Matrix) +
    scale_y_reverse() +  # Reverse the y-axis to match the matrix row order
    labs(x = "Column", y = "Row", 
         title = title) +
    theme_minimal()
  
  # Print the plot
  print(p)
}


cluster_accuracy = function(x,true_x){
  
  FDR = sum((x==1)&(true_x==0)) / sum(x==1)
  Power = sum((x==1)&(true_x==1)) / sum(true_x==1)
  Precision = mean(x==true_x)
  return(cbind(FDR,Power,Precision))
}

RowCol_from_ColMat_idx = function(n, ColMat_idx){
  # Compute the row and column indices
  row_index <- (ColMat_idx - 1) %% n + 1
  col_index <- (ColMat_idx - 1) %/% n + 1
  
  return(cbind(row_index, col_index))
}

plot_lines_by_clus = function(X,Z){
  # plot all cluster
  df <- as.data.frame(X)
  df$row_id <- 1:nrow(df)
  df$cluster = as.factor(Z)
  
  # Reshape the data frame to long format
  df_long <- melt(df, id.vars = c("row_id", "cluster"), variable.name = "Column", value.name = "Value")
  
  
  # Plot using ggplot
  ggplot(df_long, aes(x = Column, y = Value, group = row_id, color = cluster)) +
    geom_line() +
    labs(x = "Column", y = "Value", color = "Cluster") +
    theme_minimal()
}

# Jaccard index for clustering
jaccard_index <- function(pred_labels, true_labels) {
  if (length(true_labels) != length(pred_labels)) {
    stop("Length of true_labels and pred_labels must be the same")
  }
  
  # Create contingency table (confusion matrix)
  contingency <- table(true_labels, pred_labels)
  
  # Solve the assignment problem to match true and predicted labels
  assignment <- clue::solve_LSAP(contingency, maximum = TRUE)
  
  # Relabel predicted clusters based on optimal assignment
  remapped_pred_labels <- pred_labels
  unique_pred <- as.numeric(colnames(contingency))
  unique_true <- as.numeric(rownames(contingency))
  
  for (i in seq_along(unique_true)) {
    remapped_pred_labels[pred_labels == unique_pred[assignment[i]]] <- unique_true[i]
  }
  
  # Create similarity matrices for true and relabeled predicted clusters
  n <- length(true_labels)
  true_matrix <- outer(true_labels, true_labels, FUN = "==")
  pred_matrix <- outer(remapped_pred_labels, remapped_pred_labels, FUN = "==")
  
  # Convert to logical matrices (excluding diagonal)
  true_matrix <- true_matrix & lower.tri(true_matrix)
  pred_matrix <- pred_matrix & lower.tri(pred_matrix)
  
  # Compute intersection and union
  intersection <- sum(true_matrix & pred_matrix)
  union <- sum(true_matrix | pred_matrix)
  
  # Avoid division by zero
  if (union == 0) return(0)
  
  # Compute Jaccard Index
  return(intersection / union)
}

toBinary <- function(num, width) {
  if (width <= 0) return("")
  result <- ""
  for (i in seq_len(width)) {
    result <- paste0(num %% 2, result)
    num <- num %/% 2
  }
  return(result)
}

dyadicLeftMapping <- function(x) {
  if (!is.numeric(x) || any(x != floor(x)) || any(x < 1)) {
    stop("Input must be an integer or vector of integers >= 1.")
  }
  
  sapply(x, function(val) {
    # Determine the length L of the binary string (including the appended "0")
    L <- ceiling(log2(val + 1))
    # Number of strings with length less than L is 2^(L-1)-1.
    prev_count <- 2^(L - 1) - 1
    # 0-based offset within the strings of length L.
    offset <- val - prev_count - 1
    
    # Convert the offset to a binary string of width L-1.
    bin_str <- toBinary(offset, L - 1)
    
    # Append "0" to indicate a left split.
    paste0(bin_str, "0")
  })
}


plot_cluster_mean = function(X,Z_est,title = NULL,quantile = T, y_upper_lim = NULL){
  df <- as.data.frame(X)
  df$cluster <- Z_est
  
  # If your X columns are not already named in a way that reflects location, assign names:
  # For example, if you want the locations to be 1, 2, …, then:
  if(is.null(colnames(df)[-ncol(df)])) {
    colnames(df)[1:(ncol(df)-1)] <- 1:(ncol(df)-1)
  }
  
  # Reshape the data from wide to long format.
  # Each row now corresponds to one observation at a specific location.
  df_long <- df %>%
    pivot_longer(
      cols = -cluster,
      names_to = "location",
      values_to = "value"
    ) %>%
    # Remove non-numeric characters (like "X") and convert to numeric
    mutate(location = as.numeric(gsub("[^0-9]", "", location)))
  
  
  # Optional: Filter clusters with more than 10 observations
  clusters_keep <- names(which(table(df$cluster) > 10))
  df_long <- df_long %>% filter(cluster %in% clusters_keep)
  
  # Compute summary statistics (mean, lower 2.5% quantile, upper 97.5% quantile)
  summary_df <- df_long %>%
    group_by(cluster, location) %>%
    summarise(
      mean_val = mean(value),
      lower = quantile(value, probs = 0.25),
      upper = quantile(value, probs = 0.75),
      .groups = 'drop'
    )
  
  if(is.null(y_upper_lim )){
    y_upper_lim = max(summary_df$mean_val) + 0.1
    
  }
  
  # Create the ggplot
  if(quantile){
    p <- ggplot(summary_df, aes(x = location, group = factor(cluster))) +
      # Shaded area between quantile bounds
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(cluster)), 
                  alpha = 0.2, color = NA) +
      # Plot mean as a solid line
      geom_line(aes(y = mean_val, color = factor(cluster)), size = 1) +
      scale_y_continuous(limits = c(0,y_upper_lim))+
      labs(
        title = title,
        x = "Location",
        y = "Cluster Mean",
        color = "Cluster",
        fill = "Cluster"
      ) +
      theme_minimal()
  }else{

    p <- ggplot(summary_df, aes(x = location, group = factor(cluster))) +
      # Plot mean as a solid line
      geom_line(aes(y = mean_val, color = factor(cluster)), size = 1) +
      scale_y_continuous(limits = c(0,y_upper_lim))+
      labs(
        title = title,
        x = "Location",
        y = "Cluster Mean",
        color = "Cluster",
        fill = "Cluster"
      ) +
      theme_minimal()
  }
    
  
  print(p)
}

plot_cluster_mean_base = function(X,Z, title = NULL,y_upper_lim = NULL){
  if(is.null(y_upper_lim) ){
    y_upper_lim = max(X) + 0.1
  }
  clus_label = which(table(Z)>10)
  clus_label = as.numeric(names(clus_label))
  i=0
  for(k in clus_label){
    if(i==0){
      plot(apply(X[Z==k,],2,mean),type = "l",col = k+1, ylim = c(0,y_upper_lim) ,
           main = title,  # Title of the plot
           ylab = "Cluster mean",      # Label for the x-axis
           xlab = "Location")
    }else{
      lines(apply(X[Z==k,],2,mean),col = k+1)
    }
    i = i+1
  }
  legend("topright",legend = paste0(clus_label,":",table(Z)[table(Z)>10]),
         col = clus_label+1,lty = 1)
  
  
  
}

plot_all_cov = function(cor_tree){
  Z_est = apply(cor_tree$mcmc$Z[,-c(total_iter-burnin)],1,mean)
  Z_tab = table(Z_est)
  Z_show = Z_tab[Z_tab>10]
  identified_cluster = as.numeric(names(Z_show))
  
  Sigma_all_clus = array(NA, dim = c(dim(cor_tree$mcmc$Sigma_inv[[1]])))
  n_cov_mcmc = length(cor_tree$mcmc$Sigma_inv)-1
  Sigma_inv_mcmc = array(NA, dim = c(n_cov_mcmc, dim(cor_tree$mcmc$Sigma_inv[[1]])))
  for(m in 1:n_cov_mcmc){
    Sigma_inv_mcmc[m,,,] = cor_tree$mcmc$Sigma_inv[[m]]
  }
  
  for(k in identified_cluster){
    Sigma_inv_mean = apply(Sigma_inv_mcmc[,,,k],2:3,mean)
    Sigma_all_clus[,,k] = chol2inv(chol(Sigma_inv_mean))
  }
  library(gridExtra)
  plot_list <- lapply(identified_cluster, function(i)  plot_mat(Sigma_all_clus[,,i],title = paste0("Sigma",i)))
  grid.arrange(grobs = plot_list, ncol = 2)
}



# resconstruct density from tree splitting probablities -------------------


# Helper function: convert vectorized splits to list structure
vector_to_splits <- function(splits_vec) {
  # Determine the depth L such that 2^L - 1 equals length(splits_vec)
  n <- length(splits_vec)
  depth <- log2(n + 1)
  if (abs(round(depth) - depth) > .Machine$double.eps^0.5) {
    stop("Length of input vector is not of the form 2^L - 1")
  }
  depth <- as.integer(round(depth))
  
  splits_list <- list()
  start <- 1
  for (level in 1:depth) {
    n_level <- 2^(level - 1)  # number of nodes at this level
    splits_list[[level]] <- splits_vec[start:(start + n_level - 1)]
    start <- start + n_level
  }
  
  return(splits_list)
}

# Helper function: get binary representation (as a vector of 0's and 1's)
leaf_to_path <- function(leaf, depth) {
  # Convert leaf index to binary representation of specified 'depth'
  # Returns a vector of 0's and 1's, where 0 means left and 1 means right.
  path <- rep(0, depth)
  for (i in 1:depth) {
    # Check the least-significant bit, then shift right
    path[depth - i + 1] <- leaf %% 2
    leaf <- leaf %/% 2
  }
  return(path)
}

# Main function: given a list of splitting probabilities, reconstruct the density
reconstruct_density <- function(splits_vector) {
  splits <- vector_to_splits(splits_vector)
  depth <- length(splits)  # total number of splits
  n_leaves <- 2^depth      # number of leaves/intervals
  densities <- numeric(n_leaves)
  
  # Loop over each leaf
  for (leaf in 0:(n_leaves - 1)) {
    path <- leaf_to_path(leaf, depth)
    mass <- 1
    node_index <- 1  # Start at the root; in splits[[1]], there is only one element.
    
    # For each level, update the mass and move to the proper node index for the next level
    for (level in 1:depth) {
      # At level 'level', splits[[level]] is a vector of length 2^(level-1)
      # The splitting probability at the current node:
      p_split <- splits[[level]][node_index]
      
      if (path[level] == 0) {
        # Left branch: use p_split directly
        mass <- mass * p_split
        # In a binary tree stored in level-order, the left child of node_index is:
        node_index <- (node_index - 1) * 2 + 1
      } else {
        # Right branch: use (1 - p_split)
        mass <- mass * (1 - p_split)
        # Right child index:
        node_index <- (node_index - 1) * 2 + 2
      }
    }
    # Multiply by (2^depth) to convert mass into density (since each interval has length 2^(-depth))
    densities[leaf + 1] <- mass 
  }
  
  return(densities)
}



# plot functions for real data analysis -----------------------------------



heat_map_DNase_seq = function(X,Z, clamp=NULL){

  method_names = deparse(substitute(Z))
  library(pheatmap)
  library(viridis)
  counts_mat = X
  cluster_vec = Z
  # 1. Order your matrix and cluster labels
  ord          <- order(cluster_vec)
  mat_ord      <- log(counts_mat[ord, ] + 1)
  clusters_ord <- cluster_vec[ord]
  
  annotation_ord <- data.frame(Cluster = cluster_vec[ord])
  rownames(annotation_ord) <- rownames(mat_ord)
  
  # 2. Compute where to put gaps: after each cluster’s block
  #    table(...) gives counts per cluster; cumsum(...) gives the end‐row index
  gaps_after <- cumsum(table(clusters_ord))
  
  # 3. Call pheatmap with gaps_row
  pheatmap(
    mat_ord,
    color      = viridis(100, option = "magma", direction = -1),
    # annotation for rows
    annotation_row = annotation_ord,
    annotation_colors = list(
      Cluster = viridis(n = length(unique(cluster_vec)), option = "viridis")
      ),
    # NO row clustering—keeps your order
    cluster_rows = FALSE,
    # still cluster columns if you like
    cluster_cols = F,
    # insert horizontal gaps
    gaps_row     = gaps_after,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color  = NA,
    main           = paste(method_names,": log(count+1) by Cluster")
  )
  
  
  # # 1. Compute ordering of samples by cluster
  # ord <- order(cluster_vec)
  # # 2. Reorder the matrix & annotation
  # logM_ord <- log(counts_mat[ord, ] + 1)
  # annotation_ord <- data.frame(Cluster = cluster_vec[ord])
  # rownames(annotation_ord) <- rownames(logM_ord)
  # # Choose your viridis palette: "viridis", "magma", "plasma", "cividis", etc.
  # my_viridis <- viridis(n = 276, option = "magma",direction = -1)  
  # if(is.null(clamp)){
  #   clamp = max(logM_ord)
  # }
  # mat_clamp <- pmin(pmax(logM_ord, 0), clamp)  
  # pheatmap(
  #   mat_clamp,
  #   color          = my_viridis,
  #   annotation_row = annotation_ord,
  #   cluster_rows   = FALSE,
  #   cluster_cols   = F,
  #   show_rownames  = FALSE,
  #   show_colnames  = FALSE,
  #   main           = paste(method_names,": log(count+1) by Cluster"),
  #   annotation_colors = list(
  #     Cluster = viridis(n = length(unique(cluster_vec)), option = "viridis")
  #   )
  # )
  
  
}


cluster_average = function(X,Z, y_upper_lim=NULL){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Z = Z_pam
  counts_mat = X
  cluster_vec = Z
  inclusion = which(table(Z)>10)
  num_include = as.numeric(names(inclusion))
  include_id = which(Z %in% num_include)
  
  
  
  # method_names = deparse(substitute(Z))
  
  counts_mat = counts_mat[include_id,]
  cluster_vec = cluster_vec[include_id]
  
  colnames(counts_mat) = 1:ncol(counts_mat)
  df <- as.data.frame(counts_mat) %>%
    # carry over sample ID and cluster
    mutate(
      sample  = rownames(.),
      cluster = cluster_vec
    ) %>%
    # pivot to long form
    pivot_longer(
      cols      = -c(sample, cluster),
      names_to  = "bin",
      values_to = "count"
    ) %>%
    # if your bin labels are numeric strings, convert
    mutate(bin = as.numeric(bin))

  
  # 1. Compute summary statistics per cluster & bin
  df_stats <- df %>%
    group_by(cluster, bin) %>%
    summarise(
      mean_ct = mean(count),
      lo_3rd  = quantile(count, 0.25),
      hi_3rd  = quantile(count, 0.75),
      lo_95   = quantile(count, 0.025),
      hi_95   = quantile(count, 0.975),
      .groups = "drop"
    )
  
  df_stats$cluster = as.factor(df_stats$cluster)
  
  ylims <- c(0, max(df_stats$hi_95))
  if(!is.null(y_upper_lim)){
    ylims[2] = y_upper_lim
  }
  
  # 2. Get the number of subjects per cluster
  cluster_counts <- df %>%
    distinct(sample, cluster) %>%
    count(cluster) %>%
    # make a named vector: names = original cluster levels
    { setNames(paste0(.$cluster, ": n=", .$n), .$cluster) }
  
  # 2. Plot with two ribbons + mean line
  ggplot(df_stats, aes(x = bin, y = mean_ct, group = cluster)) +
    # 95% CI ribbon (light gray, semi‐transparent)
    geom_ribbon(aes(ymin = lo_95, ymax = hi_95),
                fill  = "lightgray",
                alpha = 0.6) +
    # Interquartile ribbon (darker gray)
    geom_ribbon(aes(ymin = lo_3rd, ymax = hi_3rd),
                fill  = "gray30",
                alpha = 0.4) +
    # Mean profile line (color by cluster)
    geom_line( size = 0.5) +
    facet_wrap(
      ~ cluster,
      nrow   = 1,
      scales = "fixed",
      labeller = labeller(cluster = cluster_counts)
    )  +
    scale_y_continuous(limits = ylims)+
    scale_color_brewer(palette = "Set1") +
    labs(
      x     = "Genomic bin",
      y     = "Count",
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text      = element_text(face = "bold")
    )
  
  
}
