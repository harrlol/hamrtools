## importing important libraries
library(dplyr)
install.packages('Seurat')
library(Seurat)
install.packages('patchwork')
library(patchwork)
library("Matrix")
library("readr")
library("glue")

## importing the data or count matrix (use 'ReadMtx'for bundle formate matrix or 'Read10X' for tubler formate matrix by 10X, or ReadSTARsolo for star, or manually)

Mido.data <- readMM("~/pipeline/STARsolo_results/Solo.out/Gene/filtered/matrix.mtx")
rownames(Mido.data) <- read_tsv("~/pipeline/STARsolo_results/Solo.out/Gene/filtered/features.tsv", col_names=FALSE)[, 2, drop=TRUE]
colnames(Mido.data) <- read_tsv("~/pipeline/STARsolo_results/Solo.out/Gene/filtered/barcodes.tsv", col_names=FALSE)[, 1, drop=TRUE]

#--or: Mido.data <- ReadMtx(mtx ="~/genome/matrix/matrix.mtx", cells="~/genome/matrix/barcodes.tsv", features="~/genome/matrix/features.tsv")
#--or: Mido.data <- ReadSTARsolo(data.dir ="~/genome/matrix/")


## setting up a Seurat object and displaying it
Mido <- CreateSeuratObject(counts = Mido.data, project = "Mido", min.cells = 3, min.features = 200)
Mido

# ------------------------------------------------------------------------------------------------------------------------

##### Pre-processing ######

## Quality Control and cells selection

# percent of mitrochondrial genes
Mido[["percent.mt"]] <- PercentageFeatureSet(Mido, pattern = "^MT-")
##----END OF RSCRIPT----##

exp <- <experiment dir>
  load(glue("~/Lab/{exp}/Rdata/Mido.RData"))

# Violin plot for QC metrics
VlnPlot(Mido, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(file=glue("~/Lab/{exp}/Rdata/pre-filter.png"))
VlnPlot(Mido, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# scatter plots for features
## plot1 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "percent.mt")
## plot2 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## plot1 + plot2

# choosing cells with high quality (here, has more than 200 but less than 2500 reads and less than 5% mt genes)
Mido <- subset(Mido, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# rechecking Q
VlnPlot(Mido, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png(file=glue("~/Lab/{exp}/Rdata/post-filter.png"))
VlnPlot(Mido, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

## plot1 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "percent.mt")
## plot2 <- FeatureScatter(Mido, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## plot1 + plot2

#------------------------

## Normalizing

Mido <- NormalizeData(Mido, normalization.method = "LogNormalize", scale.factor = 10000)

#------------------------

## Features selection (select genes of high variability)

# 2000 genes of highest variability are selected
Mido <- FindVariableFeatures(Mido, selection.method = "vst", nfeatures = 2000)

# these are the top 10 of them to show in the plot below
top10 <- head(VariableFeatures(Mido), 10)

# plot the variable features
plot1 <- VariableFeaturePlot(Mido)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
top10

#-----------------------

## Linear Transformation (Scaling)

all.genes <- rownames(Mido)
Mido <- ScaleData(Mido, features = all.genes)

#----------------------------------------------------------------------------------------------------------------------------

###### Linear Dimensional Reduction with PCA ####### 

## performing PCA
Mido <- RunPCA(Mido, features = VariableFeatures(object = Mido))

## examining PCA results
##print(Mido[["pca"]], dims = 1:5, nfeatures = 5) 

## visualizing PCA results in VizDimReduction(), DimPlot(), and DimHeatmap(), the heat map is particularly useful for downstream analysis
## VizDimLoadings(Mido, dims = 1:2, reduction = "pca")
## DimPlot(Mido, reduction = "pca")
## DimHeatmap(Mido, dims = 1:12, cells = 500, balanced = TRUE)

###### Determine the dataset dimentionality ########

## determining the top PCs (and how many) that represent a robust compression for the dataset.
##Mido <- JackStraw(Mido, num.replicate = 100, dims = 20)
##Mido <- ScoreJackStraw(Mido, dims = 1:20)

# Visualizing the distribution of p-values for each PC with uniform distribution (the dashed line). Elbow plot can be used as an alternative.
##JackStrawPlot(Mido, dims = 1:20)
ElbowPlot(Mido, ndims=50)
png(file=glue("~/Lab/{exp}/Rdata/elbowPlot.png"))
ElbowPlot(Mido, ndims=50)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------

####### Cells Clustering #######

## K- nearest neighbor
Mido <- FindNeighbors(Mido, dims = 1:25)

## Clustering
Mido <- FindClusters(Mido, resolution = 0.5)

## Showing the cluster ID of the first few cells in the dataset (the function Idents() calls the cluster ID)
head(Idents(Mido), 20)


## Optional: non-linear dimentionality reduction (as UMAP and tSNE) to visualize and explore the data
Mido <- RunUMAP(Mido, dims = 1:25)
DimPlot(Mido, reduction = "umap", label=TRUE)
png(file=glue("~/Lab/{exp}/Rdata/clusters.png"))
DimPlot(Mido, reduction = "umap", label=TRUE)
dev.off()

## view cell counts in each cluster
table(Mido@meta.data$seurat_clusters)
write.table(table(Mido@meta.data$seurat_clusters), 
            file=glue("~/Lab/{exp}/Description.txt"), 
            append=TRUE,
            sep="\t",
            quote=FALSE)

## Save progress to retrieve it if needed for other downstream analysis as ones not included here

saveRDS(Mido, file = glue("~/Lab/{exp}/Rdata/Mido_Seurat.rds"))

## Log top genes of each cluster
Mido.markers <- FindAllMarkers(Mido, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Mido.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(
  Mido.markers %>% group_by(cluster) %>% slice_max(n = 4, order_by = avg_log2FC), 
  file=glue("~/Lab/{exp}/Description.txt"), 
  append=TRUE, 
  sep="\t", 
  quote=FALSE)
# -----------------------------------------------------------------------------------------------------------------------------

######## Obtaining a dataframe of Cells barcodes and their corresponding cluster ID ############

CBs_Clusters_dataframe <- FetchData(Mido, vars = 'ident')
## save it as 
write.csv(CBs_Clusters_dataframe, glue("~/Lab/{exp}/Rdata/CBs_Clusters_dataframe.csv"))

## exit R environment
q()

#------------------------

