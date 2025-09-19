dropbox_url = "https://www.dropbox.com/scl/fo/usaeb9iwcomugrjulqkql/AKKKMpPTZsRlhwYK53h9D-A/processed_data/hg38?dl=0&rlkey=n841q5hun1v09z7cc47zfl9wg&subfolder_nav_tracking=1"

data_path = "/hpc/group/mastatlab/yx306/DNA_data/"
list.files(path = data_path)

source("./R/microbiome_help.R")
file_name_K_new = "NRF1.K562.DNase.counts.mat.rds"
label_K_new = "NRF1.K562.sites.chip.labels.rds"



dnase_count_matrix <- readRDS(paste0(data_path,file_name_K_new))
sites_chip_labels <- readRDS(paste0(data_path,label_K_new))
dnase_count_matrix <- dnase_count_matrix[,1:(ncol(dnase_count_matrix)/2)] + dnase_count_matrix[,(ncol(dnase_count_matrix)/2+1):ncol(dnase_count_matrix)]
pos_idx <- which(sites_chip_labels$chip_label == 1)
neg_idx <- which(sites_chip_labels$chip_label == 0)


threshold = 50
pwm_score_cutoff = 13
idx_include = which(rowSums(dnase_count_matrix) >= threshold & 
                      sites_chip_labels$pwm.score >= pwm_score_cutoff &
                      sites_chip_labels$chr %in% c("chr1"))



length(idx_include)

set.seed(2025)
all_include = idx_include
X = as.matrix(dnase_count_matrix[all_include,]); dim(X)
chip_labels = as.numeric(sites_chip_labels$chip_label[all_include])
table(chip_labels)


pos_idx = which(chip_labels ==1)
neg_idx = which(chip_labels ==0)

plot_pos_neg_profiles(X, chip_labels,
                      title = "Aggregate DNase-seq profiles around REST motifs in K562")
summary(rowSums(X))



n_clus = 5
y_upper_lim = 15
in_tree_depth = 9

pwm_score = sites_chip_labels$pwm.score[all_include]

quantile_groups <- function(x, k) {
  qs <- quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE)
  cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    labels = seq_len(k)
  )
}

initial_value = as.numeric(quantile_groups(pwm_score,k=n_clus))
table(initial_value,chip_labels)
mclust::adjustedRandIndex(initial_value,chip_labels)

library(cluster)

n_clus_bin=5
# Perform k-means clustering ------------------------------------------------




data_matrix <- scale(X)
kmeans_result <- kmeans(data_matrix, centers = n_clus_bin, nstart = 25)
Z_kmeans = kmeans_result$cluster;
plot_pos_neg_profiles(X,Z_kmeans,title="kmeans")



# PAM -----------------------------------------------------------------------



pamx <- pam(X, n_clus_bin)
Z_pam = c(pamx$clustering)


plot_pos_neg_profiles(X,Z_pam,title="PAM")

table(Z_pam,chip_labels)
mclust::adjustedRandIndex(Z_pam,chip_labels)






# try CENTIPEDE -----------------------------------------------------------

# install.packages("CENTIPEDE",
#                  repos="http://R-Forge.R-project.org")
library(CENTIPEDE)



centFit <- fitCentipede(
  Xlist = list(DNase = as.matrix(X)),
  Y     = cbind(1, pwm_score)
  
)

# posterior binding probabilities (i.e. cluster assignments)
post_probs <- centFit$PostPr


# if you want a hard 2-way cluster assignment:
cluster_assignments <- ifelse(post_probs > 0.5, 1, 0)
Z_centipede = cluster_assignments

# mclust::adjustedRandIndex(chip_labels, Z_centipede)
# run cor-tree ------------------------------------------------------------


tree_depth <- in_tree_depth # final result
print(paste("num of leaf nodes",2^(tree_depth-1))) 
cutoff_layer <- 3; print(paste("num of cor nodes",2^(cutoff_layer+1)-1))    # Layer for correlated nodes cutoff
warm_start <- 0
burnin <- warm_start + 20  # Burn-in period
total_iter <- burnin + 20;print(total_iter)  # Total iterations
c_sigma2_vec <- 10  # Hyperparameter for sigma^2 vector. 1/c is the mean
sigma_mu2 <- 0.1     # Hyperparameter for sigma_mu^2
# cov_interval <- 5 # for data_G only
cov_interval <- 3 # for data_G only
# Call the CorTree_sampler function
devtools::load_all()
set.seed(2025)
result <- CorTree_sampler(X, 
                           init_Z = Z_kmeans-1,
                           n_clus = n_clus, 
                           tree_depth, 
                           cutoff_layer, 
                           total_iter, 
                           burnin, 
                           c_sigma2_vec, 
                           sigma_mu2,
                           warm_start = warm_start,
                           cov_interval = cov_interval,
                           all_ind = F) # set sn_cor, smallest number in each cluster
# to impose correlation structure
result$elapsed
par(mfrow = c(1,1))
plot(result$mcmc$loglik[burnin:total_iter],type="l")

Z_cortree = apply(result$mcmc$Z[,-c(total_iter-burnin)],1,mean);table(Z_cortree)
Z_cortree = round(Z_cortree)
# mclust::adjustedRandIndex(chip_labels, Z_cortree)
# hist_by_group_gg(chip_data,Z_cortree, title = "Cor-tree",y_max = max_log_count)
cluster_mean_gg(X,Z_cortree, title = "Cor-tree",y_max=y_upper_lim)

saveRDS(Z_cortree,file.path(data_path,"cortree_NRF1_K_DNase_wKMeansinit.rds"))


# check convergence of the mixing proportion


n_mcmc = ncol(result$mcmc$pi)
matplot(
  t(result$mcmc$pi[,-n_mcmc]),          # now 50 rows × 5 columns → 5 “series” of length 50
  type = "l",    # plot lines
  lty  = 1,      # solid line type for all
  xlab = "Index (1:50)",
  ylab = "Value",
  main = "Convergence of mixing weights after burn-in"
)
legend(
  "topright",
  legend = paste0("Cluster ", 1:5),
  col    = 1:5,    # default color palette
  lty    = 1,
  bty    = "n"
)




# summarize with CENTIPEDE ------------------------------------------------

Z_cortree = readRDS(file.path(data_path, "cortree_NRF1_K_DNase_wKMeansinit.rds"))
Z_cortree = round(Z_cortree)

Z_cortree_show = Z_cortree
# # relabel Z_cortree for visual effect
Z_cortree_show <- ifelse(Z_cortree == 0, 1,
                         ifelse(Z_cortree == 2, 0, NA))  # NA if any other value

ari_kmeans <- mclust::adjustedRandIndex(chip_labels, Z_kmeans)
ari_pam <- mclust::adjustedRandIndex(chip_labels, Z_pam)
ari_cortree <- mclust::adjustedRandIndex(chip_labels, Z_cortree_show)
ari_initial <- mclust::adjustedRandIndex(chip_labels, initial_value)
ari_centipede <- mclust::adjustedRandIndex(chip_labels, Z_centipede)

# combine into a data.frame
ari_table <- data.frame(
  Method = c("K-means", "PAM", "CENTIPEDE", "PWM_cate","Cortree"),
  ARI    = c(ari_kmeans, ari_pam, ari_centipede, ari_initial, ari_cortree)
)

saveRDS(ari_table,file.path(data_path,"NRF1_ari_table.rds"))

knitr::kable(ari_table, digits = 3,
             caption = "Adjusted Rand Index by clustering method")

# knitr::kable(ari_table, digits = 3,format="latex",
#              caption = "Adjusted Rand Index by clustering method")

table(chip_labels, Z_cortree_show)
table(chip_labels, Z_pam)

# check histogram with Chip data
chip_data = sites_chip_labels$chip[all_include]
source("./R/microbiome_help.R")
library(patchwork)

max_log_count = 200
p0 <- hist_by_group_gg(chip_data,chip_labels, title = "chip_labels",y_max = max_log_count)
p1 <- hist_by_group_gg(chip_data,Z_kmeans, title = "Kmeans",y_max = max_log_count)
p2 <- hist_by_group_gg(chip_data,Z_pam, title = "PAM",y_max = max_log_count)
p3 <- hist_by_group_gg(chip_data,Z_centipede, title = "CENTIPEDE",y_max = max_log_count)
p4 <- hist_by_group_gg(chip_data,Z_cortree_show, title = "Cor-tree",y_max = max_log_count)

# 
# combined <- (p0+p1 )/( p2 + p3)   # or: p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# print(combined)
y_upper_lim = NULL
m0 = cluster_mean_gg(X,chip_labels, title = "chip_labels",y_max=y_upper_lim)
m1 = cluster_mean_gg(X,Z_kmeans, title = "Kmeans",y_max=y_upper_lim)
m2 = cluster_mean_gg(X,Z_pam, title = "PAM",y_max=y_upper_lim)
m3 = cluster_mean_gg(X,Z_centipede, title = "CENTIPEDE",y_max=y_upper_lim)
m4 = cluster_mean_gg(X,Z_cortree_show, title = "Cor-tree",y_max=y_upper_lim)

# combined <- (m0+m1 )/( m2 + m3)   # or: p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# print(combined)

# combined <- (p0 +p1+p2+p3 )/(m0+m1+ m2 + m3)  # or: p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
# print(combined)

combined <- wrap_plots(
  m0, m1,
  m2, m3,m4,
  p0, p1,
  p2, p3,p4,
  ncol = 5
)

print(combined)

combined <- wrap_plots(
  m0, m1,
  m2, m4,
  p0, p1,
  p2, p4,
  ncol = 4
)

print(combined)


# combine ari table -------------------------------------------------------


ari_table1 = readRDS(file.path(data_path,"NRF1_ari_table.rds"))
ari_table2 = readRDS(file.path(data_path,"REST_ari_table.rds"))
df = ari_table2
colnames(df)[2] = "REST_ARI"
df[,3] = ari_table1[,2]
colnames(df)[3] = "NRF1_ARI"

df_show = t(df)
df_show = df_show[-1,]
df_show <- apply(df_show, 2, as.numeric)
# optional: preserve row names
rownames(df_show) <- c("REST_ARI", "NRF1_ARI")
colnames(df_show) = df$Method


knitr::kable(df_show, digits = 3, "latex",
             caption = "Adjusted Rand Index by clustering method")

