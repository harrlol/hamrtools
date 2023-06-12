library(dplyr)
library(Seurat)
library(patchwork)
library("Matrix")
library("readr")
library("glue")

args=commandArgs(trailingOnly=TRUE)

# Initializing Seurat project for further visualization
mtx.data <- readMM(file.path(args[1], "matrix.mtx.gz"))
rownames(mtx.data) <- read_tsv(file.path(args[1], "features.tsv.gz"), 
                               col_names=FALSE)[, 1, drop=TRUE]
colnames(mtx.data) <- read_tsv(file.path(args[1], "barcodes.tsv.gz"), 
                               col_names=FALSE)[, 1, drop=TRUE]
mtx.seu <- CreateSeuratObject(counts = mtx.data, project = "TOBESET", min.cells = 3, min.features = 200)
mtx.seu[["percent.mt"]] <- PercentageFeatureSet(mtx.seu, pattern = "^MT-")

# Violin plot for QC metrics
png(file=file.path(args[2], "pre-filter.png"), width = 10, height = 8, units = "in", res=1200)
VlnPlot(mtx.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Scatter plots for features
png(file=file.path(args[2],"pre-featurescatter.png"), width = 10, height = 8, units = "in", res=1200)
FeatureScatter(mtx.seu, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(mtx.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Filtering cells with high quality 
# Benchmark set:  more than 200 but less than 2500 reads and less than 5% mt genes
mtx.seu <- subset(mtx.seu, subset = nFeature_RNA > args[3] & nFeature_RNA < args[4] & percent.mt < args[5])

# Rechecking the above 2 plots for comparison
png(file=file.path(args[2],"post-filter.png"), width = 10, height = 8, units = "in", res=1200)
VlnPlot(mtx.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png(file=file.path(args[2],"post-featurescatter2.png"), width = 10, height = 8, units = "in", res=1200)
FeatureScatter(mtx.seu, feature1 = "nCount_RNA", feature2 = "percent.mt") + FeatureScatter(mtx.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Normalizing seurat object
mtx.seu <- NormalizeData(mtx.seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Select 2000 genes of highest variability
mtx.seu <- FindVariableFeatures(mtx.seu, selection.method = "vst", nfeatures = 2000)

# Tp10 included below
top10 <- head(VariableFeatures(mtx.seu), 10)
write.table(top10, 
            file=file.path(args[2],"Top10VarGenes.txt"), 
            append=TRUE,
            sep="\t",
            quote=FALSE)


# Plot the variable features
plot1 <- VariableFeaturePlot(mtx.seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(file=file.path(args[2], "varfeatures.png"), width = 10, height = 8, units = "in", res=1200)
plot1 + plot2
dev.off()

## PCA
all.genes <- rownames(mtx.seu)
mtx.seu <- ScaleData(mtx.seu, features = all.genes)
mtx.seu <- RunPCA(mtx.seu, features = VariableFeatures(object = mtx.seu))

# PCA result visualization
png(file=file.path(args[2],"dimloadings.png"), width = 10, height = 8, units = "in", res=1200)
VizDimLoadings(mtx.seu, dims = 1:2, reduction = "pca")
dev.off()

png(file=file.path(args[2],"dimplot.png"), width = 10, height = 8, units = "in", res=1200)
DimPlot(mtx.seu, reduction = "pca")
dev.off()

png(file=file.path(args[2],"dimheatmap.png"), width = 10, height = 8, units = "in", res=1200)
DimHeatmap(mtx.seu, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

# Obtain PCs 
mtx.seu <- JackStraw(mtx.seu, num.replicate = 100, dims = 20)
mtx.seu <- ScoreJackStraw(mtx.seu, dims = 1:20)

# Elbow plot
png(file=file.path(args[2],"jackstrawPlot.png"), width = 10, height = 8, units = "in", res=1200)
JackStrawPlot(mtx.seu, dims = 1:20)
dev.off()

png(file=file.path(args[2],"elbowPlot.png"),  width = 10, height = 8, units = "in", res=1200)
ElbowPlot(mtx.seu, ndims=50)
dev.off()

# K-mean clustering
mtx.seu <- FindNeighbors(mtx.seu, dims = 1:25)
mtx.seu <- FindClusters(mtx.seu, resolution = 0.5)

# UMAP reduction
mtx.seu <- RunUMAP(mtx.seu, dims = 1:25)

# Save clustering plot
png(file=file.path(args[2], "clusters.png"), width = 10, height = 8, units = "in", res=1200)
DimPlot(mtx.seu, reduction = "umap", label=TRUE)
dev.off()

# Record cell counts in each cluster
write.table(table(mtx.seu@meta.data$seurat_clusters), 
            file=file.path(args[2],"cell_per_cluster.txt"), 
            append=TRUE,
            sep="\t",
            quote=FALSE)

# Save progress to retrieve it if needed for other downstream analysis as ones not included here
saveRDS(mtx.seu, file=file.path(args[2],"mtx_seu.rds"))

# Log top genes of each cluster
mtx.seu.markers <- FindAllMarkers(mtx.seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(
  mtx.seu.markers %>% group_by(cluster) %>% slice_max(n = 4, order_by = avg_log2FC), 
  file=file.path(args[2],"top_gene_per_cluster.txt"), 
  append=TRUE, 
  sep="\t", 
  quote=FALSE)

# Obtaining a dataframe of Cells barcodes and their corresponding cluster ID
CBs_Clusters_dataframe <- FetchData(mtx.seu, vars = 'ident')
write.csv(CBs_Clusters_dataframe, file=file.path(args[2], "CBs_Clusters_dataframe.csv"))
