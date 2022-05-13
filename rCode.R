library(dplyr)
library(ggplot2)
library(Seurat)

newdat <- Read10X(data.dir = "C:/Users/Data/")

dim(newdat)
summary(colSums(newdat))
View(newdat)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = newdat, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc$nUMI <- Matrix::colSums(pbmc)

slotNames(pbmc)
head(pbmc@meta.data)
View(pbmc@meta.data)
summary(pbmc@meta.data)

summary(pbmc)
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA","nUMI", "percent.mt"), ncol = 3)
dim(pbmc) #13580 x 1803
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

##########Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##############Perform linear dimensional reduction#################

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)


################Determine the 'dimensionality' of the dataset#############
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

##################Cluster the cells#########################
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
library(umap)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

#Run non-linear dimensional reduction (UMAP/tSNE)

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = pbmc, reduction = "tsne")

#############Finding differentially expressed features (cluster biomarkers)######################

# find all markers of cluster 1
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)
print(x = head(x = cluster0.markers, n = 10))

x0 <- c("ALAS2",  "AHSP",   "SLC25A37",   "MT-ND3", "HBA2", "HBB", "B2M", "HBA1", "HLA-B","SNCA")

VlnPlot(pbmc, features = x0, slot = "counts", log = TRUE,pt.size = 0)

######################8Cluster########################
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)
x1 <- c("RP11-1143G9.4",  "MNDA",   "FCN1",   "VCAN", "CSTA", "AIF1", "CTSS", "LST1", "TYROBP","FGL2")
VlnPlot(pbmc, features = x1, slot = "counts", log = TRUE,pt.size = 0)

cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)
x2 <- c("CD8B",  "GZMK",   "CMC1",   "CD8A", "CD2", "KLRD1", "CST7", "CD3E", "IL32","GZMA")
VlnPlot(pbmc, features = x2, slot = "counts", log = TRUE,pt.size = 0)


cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 10)
x3 <- c("IL7R",  "LEF1",   "PRKCQ-AS1",   "MAL", "CCR7", "TCF7", "LTB", "LDHB", "NOSIP","CD27")
VlnPlot(pbmc, features = x3, slot = "counts", log = TRUE,pt.size = 0)

cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
x4 <- c("FGFBP2",  "ZNF683",   "GZMB",   "PRF1", "HOPX", "GZMH", "PATL2", "GNLY", "CST7","CTSW")
VlnPlot(pbmc, features = x4, slot = "counts", log = TRUE,pt.size = 0)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 10)
x5 <- c("TUBB1",  "ITGA2B",   "C6orf25", "SH3BGRL2", "PRKAR2B", "C19orf33", "CMTM5", "GP9", "SPARC","SDPR")
VlnPlot(pbmc, features = x5, slot = "counts", log = TRUE,pt.size = 0)

cluster6.markers <- FindMarkers(pbmc, ident.1 =6, min.pct = 0.25)
head(cluster6.markers, n = 10)
x6 <- c("KLRC1",  "SH2D1B",   "KLRF1",   "KLRB1", "TMIGD2", "IL2RB", "TNFRSF18", "TXK", "PRSS23","PRF1")
VlnPlot(pbmc, features = x6, slot = "counts", log = TRUE,pt.size = 0)

cluster7.markers <- FindMarkers(pbmc, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 10)
x7 <- c("VPREB3",  "MS4A1",   "CD79A",   "BANK1", "LINC00926", "RALGPS2", "CD22", "SPIB", "FCER2","HLA-DOB")
VlnPlot(pbmc, features = x7, slot = "counts", log = TRUE,pt.size = 0)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)


write.csv(pbmc.markers,"D:/Lung/results/DEGLung.csv", row.names = FALSE)

#DoHeatmap(pbmc, features = top5$gene) + NoLegend()
dim(pbmc)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene, label = F) + NoLegend()
DoHeatmap(pbmc, features = top10$gene, label = F)+ scale_fill_gradientn(colors = c("blue", "white", "red"))
#DoHeatmap(pbmc)+ scale_fill_gradientn(colors = c("blue", "white", "red"))
################VolcanoPlot################################
res <- read.csv("D:/Lung/results/DEGLung.csv")
head(res)
dim(res)


with(res, plot(avg_logFC, -log10(p_val_adj), pch=20, main="Volcano plot", cex.axis=1.4))
with(subset(res, p_val_adj<.01 & avg_logFC >2), points(avg_logFC, -log10(p_val_adj), pch=20, col="red", cex.axis=1.4))
with(subset(res, p_val_adj<.01 & avg_logFC <(-2)), points(avg_logFC, -log10(p_val_adj), pch=20, col="green", cex.axis=1.4))


DEg <- subset(res, p_val_adj <.01 & abs(avg_logFC)>=2)

dim(DEgfdr2)

Down <- subset(DEg, p_val_adj <.01 & (avg_logFC<(-2)))
Up <- subset(DEg, p_val_adj <.01 & avg_logFC > 2)

dim(Down)
dim(Up)

write.csv(DEg, file = "D:/Lung/results/DEg.csv")
write.csv(Down, file = "D:/Lung/results/Down.csv")
write.csv(Up, file = "D:/Lung/results/Up.csv")


