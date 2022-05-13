library(monocle)
library(Rtsne)
library(dplyr)
library(igraph)
library(ggrepel)
library(Biobase)
library(BiocGenerics)
library(cellrangerRkit)
library(devtools)
library(DDRTree)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(clues)


my_dir <- "E:/Data/Lung_data"


# load data
gbm <- load_cellranger_matrix(my_dir)
dim(exprs(gbm))

my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = fData(gbm)),
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

# rename gene symbol column
my_feat <- fData(gbm)
names(my_feat) <- c('id', 'gene_short_name')

# no warning
my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = my_feat),
                         lowerDetectionLimit = 0.8,
                         expressionFamily = negbinomial.size())

slotNames(my_cds)

# normalisation and variance estimation steps
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

##################Clustering cells without marker genes########################################################
disp_table <- dispersionTable(my_cds)
head(disp_table)



#We will select genes, which have a mean expression >= 0.1, to use in the clustering step.
table(disp_table$mean_expression>=0.1)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)

plot_ordering_genes(my_cds)

# if num_dim is not set, the default is set to 50
# I'll use 5 here based on plot_pc_variance_explained
#?reduceDimension
my_cds <- reduceDimension(my_cds, max_components = 2, num_dim =10,
                          reduction_method = 'tSNE', verbose = TRUE)


# perform unsupervised clustering requesting 20 clusters
my_cds <- clusterCells(my_cds, num_clusters = 12)

#without specifying the number of clusters
my_cds <- clusterCells(my_cds)

# cluster information is in the Cluster column
head(pData(my_cds))
# store cluster info for comparing later
my_cluster_dim_5 <- pData(my_cds)$Cluster

plot_cell_clusters(my_cds)

#I'll perform the differential expression analysis on the subset of genes where the mean expression (across cells) was >= 0.1 (see section "Clustering cells without marker genes").

length(unsup_clustering_genes$gene_id)
de_cluster_one <- differentialGeneTest(my_cds[unsup_clustering_genes$gene_id,],
                                       fullModelFormulaStr = '~Cluster',
                                       cores = 8)
dim(de_cluster_one) 
de_cluster_one %>% arrange(qval) %>% head()

my_ordering_genes <- row.names(de_cluster_one)[order(de_cluster_one$qval<0.05)][1:1000]
my_cds_subset <- setOrderingFilter(my_cds, ordering_genes = my_ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')

my_cds_subset <- orderCells(my_cds_subset)
plot_cell_trajectory(my_cds_subset, color_by = "Cluster")

head(pData(my_cds_subset))

my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 8)

my_pseudotime_de %>% arrange(qval) %>% head()

# save the top 6 genes
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(id) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$id

plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])

# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval<0.05) %>% head(1000) %>% select(id) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster$id

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                 num_clusters = 7,
                                                 cores = 8,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

# hierarchical clustering was used to cluster the genes
# we can cut the dendrogram to form the same 3 clusters as plot_pseudotime_heatmap
my_cluster <- cutree(my_pseudotime_cluster$tree_row, 7)
my_cluster

# genes in cluster 1
gene <- my_pseudotime_de[names(my_cluster[my_cluster ==1]),"gene_short_name"]

clust_DE <- write.csv(gene, file = "C:/Users/HX-JTA/Desktop/fig/clustr1.csv")

#Analyzing Branches in Single-Cell Trajectories
plot_cell_trajectory(my_cds_subset, color_by = "Cluster")

