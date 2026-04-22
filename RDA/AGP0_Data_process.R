
rm(list = ls())
Tree_path = normalizePath(".")
ag_data = file.path(Tree_path, 'data')
# BiocManager::install("phyloseq")


library(data.table)
library(phyloseq)
library(ape)
library(dplyr)

AG_tree = phyloseq::read_tree_greengenes(file.path(ag_data,"97_otus.tree"))        #-import the phylogenetric tree
otu = fread(file.path(ag_data,'ag_fecal_from_biom.txt'), header = TRUE)  #-import the otu table


otu_matrix = do.call(rbind, otu)
otu_matrix = apply(otu_matrix,1, as.numeric)
row.names(otu_matrix) = otu_matrix[, "OTUID"]
otu_neat = otu_matrix[, -1]
rm(otu_matrix)


#################################################################################################
####----choose subset of samples

ag_fecal = fread(file.path(ag_data,"ag_fecal.txt"))
# ag_fecal = ag_fecal[ag_fecal$SEX == "male" | ag_fecal$SEX == "female"] # JM didn't use this

# ag_fecal = sample_n(ag_fecal, 150)
# otu_neat = otu_neat[, intersect(colnames(otu_neat), ag_fecal$SampleID) ]

#################################################################################################
####----choose subset of otus based on some criterion

IBD_diagonosed = which(ag_fecal$IBD %in% c(
  "Diagnosed by a medical professional (doctor, physician assistant)"
));length(IBD_diagonosed)



Diabetes_diagonosed = which(ag_fecal$DIABETES %in% c(
  "Diagnosed by a medical professional (doctor, physician assistant)"
));length(Diabetes_diagonosed)

KIDNEY_DISEASE_diagonosed = which(ag_fecal$KIDNEY_DISEASE %in% c(
  "Diagnosed by a medical professional (doctor, physician assistant)"
));length(KIDNEY_DISEASE_diagonosed)



chosen_sample = union(union(IBD_diagonosed, Diabetes_diagonosed),KIDNEY_DISEASE_diagonosed);
length(chosen_sample)

otu_subgroup = otu_neat[, intersect(colnames(otu_neat), ag_fecal$SampleID[ IBD_diagonosed ])]

otu_subgroup = otu_neat[, colnames(otu_neat)]
dim(otu_subgroup)
otu_rowsum = apply(otu_subgroup, 1, sum)
otu_colsum = apply(otu_subgroup, 2, sum)

# JM's setting:
num_otu_top = 75  ####----number of top otus to be used in the analysis
num_min_sample = 500

# # Try with more subjects
# num_otu_top = 64  ####----number of top otus to be used in the analysis
# num_min_sample = 100

cutoff = sort(otu_rowsum, decreasing = TRUE)[num_otu_top]

otu_top = otu_subgroup[otu_rowsum >= cutoff, ]
otu_top = otu_top[, apply(otu_top, 2, sum) >= num_min_sample]
dim(otu_top)

# ag_fecal = ag_fecal[ag_fecal$SampleID %in% colnames(otu_top), ]

# YX: check heatmap of otu_top ---------
library(reshape2)
library(ggplot2)

# Create the sample rectangular matrix (as above)
X = t(otu_top)
summary(rowSums(X))
summary(colSums(X))
data_matrix <- log(X+1)
data_matrix = data_matrix[order(rowSums(data_matrix)),]
rownames(data_matrix) <- paste("Row", 1:nrow(X), sep = "")
colnames(data_matrix) <- paste("Col", 1:ncol(X), sep = "")


# Reshape the matrix into long format
df <- melt(data_matrix)
colnames(df) <- c("Row", "Column", "Value")

# Plot the heatmap using ggplot2
ggplot(df, aes(x = Column, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  scale_x_discrete(labels = c(1:221)) +  # Set custom x-axis labels
  labs(title = "Heatmap of log(X+1). X = Count of OTUs",
       x = "OTU ID")+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),   # Hide row names
    axis.ticks.y = element_blank()   # Hide row ticks
  )
saveRDS(X, file = "./data/OTU_Diabetes.rds")


# knn ---------------------------------------------------------------------

data_matrix <- scale(X)
kmeans_result <- kmeans(data_matrix, centers = 5, nstart = 25)
Z = kmeans_result$cluster; Z_kmeans = Z
table(Z)
clus_label = which(table(Z)>10)
par(mfrow = c(1,1))
i=0
for(k in clus_label){
  if(i==0){
    plot(apply(X[Z==k,],2,mean),type = "l",col = k, ylim = c(0,max(rowMeans(X))),
         main = "k-means",  # Title of the plot
         ylab = "Cluster mean",      # Label for the x-axis
         xlab = "Location")
  }else{
    lines(apply(X[Z==k,],2,mean),col = k)
  }
  i = i+1
}
legend("topright",legend = clus_label,col = clus_label,lty = 1)

Z_clus = table(Z);max(Z_clus)/sum(Z_clus)

#################################################################################################
####---create a trimmed otu tree with only the otus selected above

diff = setdiff(AG_tree$tip.label, row.names(otu_top))
tree = drop.tip(AG_tree, tip=diff, trim.internal = TRUE, subtree = FALSE, root.edge = 0)

#################################################################################################
####----remove the ununsed variables

rm(otu, otu_neat, AG_tree, cutoff, diff, num_otu_top, otu_rowsum, otu_colsum)

#################################################################################################
####----save the data
saveRDS(list(otu_top = otu_top, tree = tree, ag_fecal = ag_fecal), file="./data/Cluster_AG_subsample.rds")
# save(otu_top, tree, ag_fecal, file="./data/Cluster_AG_subsample_IBD.RData")
