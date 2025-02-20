##### Preprocessing and QC Script
############################################## Loading in the required packages############################################## 
suppressMessages({
  library(dior)
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
  library(scater)
  library(scran)
  library(patchwork)
  library(dplyr)
  library(Seurat)
  library(SoupX)
  library(glmGamPoi)
  library(lattice)
  library(DoubletFinder)
  library(parallel) # detectCores()
  library(sceasy)
  library(anndata)
  library(cowplot)
  library(reticulate)
  library(harmony)
  library(SeuratDisk)
  library(anndata)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
set.seed(42)

############################################## Import data and prepare paths for QC outputs ############################################## 
# Loading in the .h5 outputs from CellRanger 6.1.2
cellranger_outputs_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Cellranger_Output_Files/'

SD1_counts <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD1_filtered_feature_bc_matrix.h5"))
SD1_raw_matrix  <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD1_raw_feature_bc_matrix.h5"), use.names = T)
SD2_counts <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD2_filtered_feature_bc_matrix.h5"))
SD2_raw_matrix  <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD2_raw_feature_bc_matrix.h5") ,use.names = T)
SD3_counts <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD3_filtered_feature_bc_matrix.h5"))
SD3_raw_matrix  <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD3_raw_feature_bc_matrix.h5"), use.names = T)
SD4_counts <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD4_filtered_feature_bc_matrix.h5"))
SD4_raw_matrix  <- Read10X_h5(filename = file.path(cellranger_outputs_dir, "SD4_raw_feature_bc_matrix.h5"), use.names = T)

# Define counts for each dataset and save outputs of SoupX
datasets <- list(SD1 = list(counts = SD1_counts, raw_matrix = SD1_raw_matrix, out_folder = "soupXOutFolder_SD1/"),
                 SD2 = list(counts = SD2_counts, raw_matrix = SD2_raw_matrix, out_folder = "soupXOutFolder_SD2/"),
                 SD3 = list(counts = SD3_counts, raw_matrix = SD3_raw_matrix, out_folder = "soupXOutFolder_SD3/"),
                 SD4 = list(counts = SD4_counts, raw_matrix = SD4_raw_matrix, out_folder = "soupXOutFolder_SD4/"))
dataset_names <- names(datasets)

# Define all output directories
soupx_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/Correct_PBMC_Analysis/QC_Figures/SoupX_Figures/"
doublet_finder_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/Correct_PBMC_Analysis/QC_Figures/Doublet_Finder_Figures/'
data_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset'

# Set directory for all subsequent outputs of SoupX
setwd(soupx_dir)

############################################## Visualize data prior to running QC pipeline ############################################## 
for (i in seq_along(datasets)) {
  dataset_name <- dataset_names[i]
  dataset <- datasets[[dataset_name]]
  
  # Create Seurat object
  srat <- CreateSeuratObject(counts = dataset$counts)
  
  # Process and visualize data prior to SoupX
  srat <- NormalizeData(srat, normalization.method = 'LogNormalize')
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat <- ScaleData(srat, vars.to.regress = c('nCount_RNA'))
  srat <- RunPCA(srat, npcs = 50, verbose = FALSE)
  srat <- RunUMAP(srat, dims = 1:50, verbose = FALSE)
  srat <- RunTSNE(srat, dims = 1:50, verbose = FALSE)
  srat <- FindNeighbors(srat, dims = 1:50, verbose = FALSE)
  srat <- FindClusters(srat, verbose = TRUE)
  srat[['percent.mt']] <- PercentageFeatureSet(srat, pattern = '^MT-')
  
  # Save UMAP and tSNE plots
  pdf(file = paste0(dataset_name, "_UMAP_tSNE_before_SoupX.pdf"))
  print(DimPlot(srat, reduction = 'tsne', group.by = 'seurat_clusters', pt.size = 0.1, label = TRUE) + NoLegend())
  print(DimPlot(srat, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.1, label = TRUE) + NoLegend())
  print(DimPlot(srat, reduction = 'tsne', group.by = 'seurat_clusters', pt.size = 0.1, label = TRUE))
  print(DimPlot(srat, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.1, label = TRUE))
  print(FeaturePlot(srat, features = c('CD8A', 'CD4', 'LYZ', 'MS4A1', 'GNLY'), order = TRUE))  # Marker genes
  dev.off()
  
  # Save QC plots
  pdf(file = paste0(dataset_name, "_QC_plot_before_SoupX.pdf"))
  print(VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, group.by = 'orig.ident', ncol = 3))
  print(densityplot((srat@meta.data$nFeature_RNA), xlim = c(0, 2000), na.rm = TRUE, plot.points = FALSE, main = 'nFeature_RNA_Dist'))
  print(densityplot((srat@meta.data$nCount_RNA), xlim = c(0, 6000), na.rm = TRUE, plot.points = FALSE, main = 'nCount_RNA_Dist'))
  print(densityplot((srat@meta.data$percent.mt), xlim = c(0, 50), na.rm = TRUE, plot.points = FALSE, main = 'percent.mt_Dist'))
  dev.off()
}

############################################## Performing ambient RNA removal using SoupX ############################################## 
for (i in seq_along(datasets)) {
  dataset_name <- dataset_names[i]
  dataset <- datasets[[dataset_name]]
  
  # Run SoupX algorithm on each sample
  srat <- CreateSeuratObject(counts = dataset$counts)
  srat <- NormalizeData(srat, normalization.method = 'LogNormalize')
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
  srat <- ScaleData(srat, vars.to.regress = c('nCount_RNA'))
  srat <- RunPCA(srat, npcs = 50, verbose = FALSE)
  srat <- RunUMAP(srat, dims = 1:50, verbose = FALSE)
  srat <- FindNeighbors(srat, dims = 1:50, verbose = FALSE)
  srat <- FindClusters(srat, verbose = TRUE)
  
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  
  soup.channel <- SoupChannel(dataset$raw_matrix, dataset$counts)
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  soup.channel <- autoEstCont(soup.channel)
  adj.matrix <- adjustCounts(soup.channel, roundToInt = TRUE)
  
  # Save adjusted counts
  DropletUtils:::write10xCounts(dataset$data_dir, adj.matrix)
  sratSoupx <- CreateSeuratObject(counts = adj.matrix)
  
  srat@meta.data$Ambient_RNA_Removal <- "Before_SoupX"
  sratSoupx@meta.data$Ambient_RNA_Removal <- "After_SoupX"
  
  # Ensure unique cell names before merging
  srat <- RenameCells(srat, add.cell.id = "noAmb")
  sratSoupx <- RenameCells(sratSoupx, add.cell.id = "Amb")
  
  # Merge objects
  Merge <- merge(x = srat, y = sratSoupx)
  
  # Ensure factor levels are set correctly
  Merge@meta.data$Ambient_RNA_Removal <- factor(Merge@meta.data$Ambient_RNA_Removal, levels = c("Before_SoupX", "After_SoupX"))
  
  # Compute mitochondrial percentage
  Merge[['percent.mt']] <- PercentageFeatureSet(Merge, pattern = "^MT-")
  
  # Set identities
  Idents(Merge) <- Merge@meta.data$Ambient_RNA_Removal
  
  # Generate violin plots
  pdf(file = paste0("Counts_Comparison_", dataset_name, ".pdf"))
  print(VlnPlot(Merge, "nCount_RNA", y.max = 10000, group.by = "Ambient_RNA_Removal"))
  print(VlnPlot(Merge, "nFeature_RNA", y.max = 3000, group.by = "Ambient_RNA_Removal"))
  print(VlnPlot(Merge, "percent.mt", y.max = 10, group.by = "Ambient_RNA_Removal"))
  dev.off()
  
  # Visualize the output of SoupX algorithm
  sratSoupx <- NormalizeData(sratSoupx, normalization.method = "LogNormalize")
  sratSoupx <- FindVariableFeatures(sratSoupx, selection.method = "vst", nfeatures = 2000)
  sratSoupx <- ScaleData(sratSoupx, vars.to.regress = c("nCount_RNA"))
  sratSoupx <- RunPCA(sratSoupx, npcs = 50, verbose = FALSE)
  sratSoupx <- RunUMAP(sratSoupx, dims = 1:50, verbose = FALSE)
  sratSoupx <- RunTSNE(sratSoupx, dims = 1:50, verbose = FALSE)
  sratSoupx <- FindNeighbors(sratSoupx, dims = 1:50, verbose = FALSE)
  sratSoupx <- FindClusters(sratSoupx, verbose = TRUE)
  sratSoupx[['percent.mt']] <- PercentageFeatureSet(sratSoupx, pattern = "^MT-")
  
  pdf(file = paste0("QC_after_SoupX_", dataset_name, ".pdf"))
  print(VlnPlot(sratSoupx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, group.by = "orig.ident", ncol = 3))
  print(FeatureScatter(sratSoupx, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident") + ggtitle("nCounts_by_percent_mt"))
  print(densityplot((sratSoupx@meta.data$nFeature_RNA), xlim = c(0, 2000), na.rm = TRUE, plot.points = FALSE, main = "nFeature_RNA_Dist"))
  print(densityplot((sratSoupx@meta.data$nCount_RNA), xlim = c(0, 6000), na.rm = TRUE, plot.points = FALSE, main = "nCount_RNA_Dist"))
  print(densityplot((sratSoupx@meta.data$percent.mt), xlim = c(0, 50), na.rm = TRUE, plot.points = FALSE, main = "percent.mt_Dist"))
  print(DimPlot(sratSoupx, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.1, label = TRUE))
  print(DimPlot(sratSoupx, reduction = "tsne", group.by = "seurat_clusters", pt.size = 0.1, label = TRUE))
  print(FeaturePlot(sratSoupx, features = c("CD8A", "CD4", "LYZ", "MS4A1", "GNLY"), order = TRUE))  # Marker genes
  dev.off()
}

############################################## Saving the outputs of SoupX ############################################## 
for (i in seq_along(datasets)) {
  dataset_name <- dataset_names[i]
  output_folder <- datasets[[dataset_name]]$out_folder
  
  # Read saved data
  adj_data <- Read10X(output_folder)
  srat <- CreateSeuratObject(counts = adj_data, project = dataset_name)
  
  saveRDS(srat, file = file.path(data_dir, paste0(dataset_name, "_srat.rds")))
}

############################################## Performing doublet removal using DoubletFinder ############################################## 

# Read in the data
SD1_srat <- readRDS(file = file.path(data_dir, "SD1_srat.rds"))
SD2_srat <- readRDS(file = file.path(data_dir, "SD2_srat.rds"))
SD3_srat <- readRDS(file = file.path(data_dir, "SD3_srat.rds"))
SD4_srat <- readRDS(file = file.path(data_dir, "SD4_srat.rds"))


# Define datasets and their sample names
datasets <- list(SD1_srat = "HC-unstim",
                 SD2_srat = "HC-stim",
                 SD3_srat = "SA-unstim",
                 SD4_srat = "SA-stim")

# Add orig.ident to meta.data for each sample
for (dataset_name in names(datasets)) {
  srat <- get(dataset_name)  # Get the Seurat object
  orig_ident_value <- datasets[[dataset_name]]  # Retrieve the orig.ident value
  
  # Assign orig.ident value
  srat@meta.data[1] <- orig_ident_value
  metadata <- srat@meta.data
  colnames(metadata)[1] <- "orig.ident"
  srat@meta.data <- metadata
  
  # Update the object in the environment
  assign(dataset_name, srat)
}

# Merge all samples into a single Seurat object
all.samples <- merge(SD1_srat, y = list(SD2_srat, SD3_srat, SD4_srat), 
                     add.cell.ids = c("SD1", "SD2", "SD3", "SD4"))

# Assign and validate identities
unique(sapply(X = strsplit(colnames(all.samples), split = "_"), FUN = "[", 1))
Idents(all.samples) <- gsub("_.*", "", colnames(all.samples))
Idents(all.samples) <- 'orig.ident'
table(all.samples@meta.data[["orig.ident"]])

# Plot QC metrics prior to DoubletFinder
all.samples$mitoPercent <- PercentageFeatureSet(all.samples, pattern = '^MT-')

# Set working directory for DoubletFinder outputs
setwd(doublet_finder_dir)

pdf('QC_plots_before_doubletfinder.pdf')
print(VlnPlot(all.samples, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3, pt.size = 0))
plot1 <- FeatureScatter(all.samples, feature1 = "nCount_RNA", feature2 = "mitoPercent")
plot2 <- FeatureScatter(all.samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(densityplot((all.samples@meta.data$nFeature_RNA), xlim = c(0, 2000), na.rm = TRUE, plot.points = FALSE, main = 'nFeature_RNA_Dist'))
print(densityplot((all.samples@meta.data$nCount_RNA), xlim = c(0, 6000), na.rm = TRUE, plot.points = FALSE, main = 'nCount_RNA_Dist'))
print(densityplot((all.samples@meta.data$mitoPercent), xlim = c(0, 50), na.rm = TRUE, plot.points = FALSE, main = 'percent.mt_Dist'))
print(plot1 + plot2)
dev.off()

# Process the merged object
all.samples <- NormalizeData(all.samples, normalization.method = "LogNormalize", scale.factor = 10000)
all.samples <- FindVariableFeatures(all.samples, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(all.samples)
all.samples <- ScaleData(all.samples, features = all.genes)
all.samples <- RunPCA(all.samples, npcs = 50, verbose = FALSE)
all.samples <- FindNeighbors(all.samples, dims = 1:40, reduction = 'pca', verbose = FALSE)
all.samples <- RunUMAP(all.samples, dims = 1:40)
all.samples <- FindClusters(all.samples, verbose = FALSE)

# Save UMAP plot
pdf('merged_data_before_doublet_finder.pdf')
print(DimPlot(all.samples, reduction = 'umap', group.by = 'orig.ident'))
dev.off()

# Save the combined Seurat object
saveRDS(all.samples, file = file.path(data_dir, 'Combined_Shrimp_Allergy_before_filtering_and_DoubletFinder.rds'))

# Read back the saved file if needed
all.samples <- readRDS(file = file.path(data_dir, 'Combined_Shrimp_Allergy_before_filtering_and_DoubletFinder.rds'))

# Run the DoubletFinder algorithm to remove doublets from all four samples
all.samples.split <- SplitObject(all.samples, split.by = "orig.ident") 
for (i in 1:length(all.samples.split)) {
  print(paste0("orig.ident",i))
  PBMC_sample    <- NormalizeData(all.samples.split[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  PBMC_sample    <- FindVariableFeatures(PBMC_sample, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(PBMC_sample)
  PBMC_sample    <- ScaleData(PBMC_sample, features = all.genes)
  PBMC_sample    <- RunPCA(PBMC_sample, npcs = 50, verbose = F)
  
  stdv <- PBMC_sample[["pca"]]@stdev
  sum.stdv <- sum(PBMC_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  PBMC_sample <- RunUMAP(PBMC_sample, dims = 1:min.pc)
  PBMC_sample <- FindNeighbors(object = PBMC_sample, dims = 1:min.pc)              
  PBMC_sample <- FindClusters(object = PBMC_sample)
  
  sweep.list <- paramSweep(PBMC_sample, PCs = 1:min.pc, sct = F)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  annotations <- PBMC_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(0.035*nrow(PBMC_sample@meta.data)) ## Assuming 3.5% according to 10X recommendations
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  PBMC_sample <- doubletFinder(seu = PBMC_sample, 
                                  PCs = 1:min.pc, 
                                  pK = optimal.pk,
                                  nExp = nExp.poi.adj,
                                  sct = F)
  metadata <- PBMC_sample@meta.data
  colnames(metadata)[8] <- "doublet"
  PBMC_sample@meta.data <- metadata 
  
  # Subset and save
  table(PBMC_sample@meta.data$doublet)
  PBMC_sample.singlets <- subset(PBMC_sample, doublet == "Singlet")
  all.samples.split[[i]] <- PBMC_sample.singlets
  remove(PBMC_sample.singlets)
}

# Merge samples into a new list containing only singlets
table(all.samples.split$`HC-stim`@meta.data$doublet)

PBMC_sample.singlets <- merge(
  x = all.samples.split[[1]],
  y = all.samples.split[-1],
  project = "shrimp_pbmc"
)

# Save and reload the singlets data
saveRDS(file = file.path(data_dir, PBMC_sample.singlets, 'Combined_after_ambient_and_doublet_removal.rds'))
PBMC_sample.singlets <- readRDS('Combined_after_ambient_and_doublet_removal.rds')
PBMC_sample.singlets <- readRDS(data_dir, 'Combined_after_ambient_and_doublet_removal.rds')

# Cluster data to compare before and after doublet removal
print(dim(PBMC_sample.singlets))  # Check cell count in object
print(dim(PBMC_sample.singlets@meta.data))  # Check number of rows in metadata
print(dim(PBMC_sample.singlets@assays$RNA@data))  # Check number of cells in expression data

PBMC_sample.singlets <- PBMC_sample.singlets %>%
  NormalizeData(normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(vars.to.regress = 'nCount_RNA') %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:50, verbose = FALSE) %>%
  RunTSNE(dims = 1:50, verbose = FALSE) %>%
  FindNeighbors(dims = 1:50, verbose = FALSE) %>%
  FindClusters(verbose = TRUE)

# Save clustering plots
pdf('Clustering_after_doubletfinder.pdf')
DimPlot(PBMC_sample.singlets, reduction = 'umap', group.by = 'seurat_clusters')
DimPlot(PBMC_sample.singlets, reduction = 'umap', group.by = 'orig.ident')
DimPlot(PBMC_sample.singlets, reduction = 'tsne', group.by = 'seurat_clusters')
dev.off()

############################################## Filtering low quality cells ############################################## 
PBMC_sample.singlets$mitoPercent <- PercentageFeatureSet(PBMC_sample.singlets, pattern = '^MT-')

# Set `orig.ident` and assign colors
PBMC_sample.singlets@meta.data$orig.ident <- factor(
  PBMC_sample.singlets@meta.data$orig.ident, 
  levels = c("HC-unstim", "HC-stim", "SA-unstim", "SA-stim")
)
my_colors <- c('#9FBDD3', '#5A6EDC', '#F6CC93', '#F09924')
PBMC_sample.singlets@meta.data$group.colors <- my_colors[as.integer(PBMC_sample.singlets@meta.data$orig.ident)]

# Generate QC plots
pdf('QC_plots_after_DoubletFinder.pdf')
VlnPlot(PBMC_sample.singlets, 
        features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), 
        ncol = 3, group.by = 'orig.ident', pt.size = 0, y.max = c(4000, 10000, 10))
for (feature in c("nFeature_RNA", "nCount_RNA", "mitoPercent")) {
  VlnPlot(PBMC_sample.singlets, 
          features = feature, 
          ncol = 1, group.by = 'orig.ident', 
          pt.size = 0, 
          y.max = ifelse(feature == "nFeature_RNA", 4000, ifelse(feature == "nCount_RNA", 10000, 50)), 
          cols = my_colors)
}
plot1 <- FeatureScatter(PBMC_sample.singlets, group.by = 'orig.ident', feature1 = "nCount_RNA", feature2 = "mitoPercent", cols = my_colors)
plot2 <- FeatureScatter(PBMC_sample.singlets, group.by = 'orig.ident', feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = my_colors)
densityplot((PBMC_sample.singlets@meta.data$nFeature_RNA), xlim = c(0, 2000), na.rm = TRUE, plot.points = FALSE, main = 'nFeature_RNA_Dist')
densityplot((PBMC_sample.singlets@meta.data$nCount_RNA), xlim = c(0, 6000), na.rm = TRUE, plot.points = FALSE, main = 'nCount_RNA_Dist')
densityplot((PBMC_sample.singlets@meta.data$mitoPercent), xlim = c(0, 50), na.rm = TRUE, plot.points = FALSE, main = 'percent.mt_Dist')
print(plot1)
print(plot2)
dev.off()

# Filter low-quality cells
PBMC_filtered <- subset(
  PBMC_sample.singlets, 
  subset = nCount_RNA > 500 & nCount_RNA < 20000 &
    nFeature_RNA > 200 & nFeature_RNA < 4000 &
    mitoPercent < 8
)

# Save Filtered Data
saveRDS(PBMC_filtered, 'Shrimp_allergy_final_after_QC.rds')
PBMC_filtered <- readRDS('Shrimp_allergy_final_after_QC.rds')

# Save in h5ad format
sceasy::convertFormat(
  PBMC_filtered, 
  from = "seurat", to = "anndata", 
  outFile = 'PBMC_filtered.h5ad'
)

############################################## Apply harmony batch correction for the different sample groups ############################################## 
PBMC_filtered <- PBMC_filtered %>%
  NormalizeData(normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(PBMC_filtered)) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:50, verbose = FALSE) %>%
  RunTSNE(dims = 1:50, verbose = FALSE) %>%
  FindNeighbors(dims = 1:50, verbose = FALSE) %>%
  FindClusters(verbose = TRUE)

# Assign colors for downstream analysis
PBMC_filtered@meta.data$orig.ident <- factor(
  PBMC_filtered@meta.data$orig.ident, 
  levels = c("HC-unstim", "HC-stim", "SA-unstim", "SA-stim")
)
PBMC_filtered@meta.data$group.colors <- my_colors[as.integer(PBMC_filtered@meta.data$orig.ident)]

# Save and plot clustering results
saveRDS(PBMC_filtered, 'Before_batch_correction_data.rds')

pdf('Clustering_after_QC.pdf')
DimPlot(PBMC_filtered, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.1, label = TRUE)
DimPlot(PBMC_filtered, reduction = 'umap', group.by = 'orig.ident', cols = my_colors, pt.size = 0.1)
DimPlot(PBMC_filtered, reduction = 'tsne', group.by = 'seurat_clusters', pt.size = 0.1)
dev.off()

# Define a helper function
create_feature_dotplots <- function(data, features, group.by, filename, point.size = 0.7) {
  pdf(file = file.path(ann_dir, filename))
  FeaturePlot(data, features = features, pt.size = point.size, order = TRUE)
  plot_grid(
    DotPlot(data, features = features, group.by = group.by) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) + 
      guides(fill = FALSE)
  )
  dev.off()
}

# Marker Gene Analysis (Broad)
create_feature_dotplots(
  PBMC_filtered, 
  features = c('CD4', 'IL7R', 'LTB', 'MAL', 'TPT1', 'LDHB', 'TRAC', 'TMSB10', 'CD3G', 'CD3D'), 
  group.by = 'Final_annotation_broad', 
  filename = 'All_markers_dotplot.pdf'
)

# Define marker groups
marker_groups <- list(
  CD4_T = c('CSF2', 'CCR10', 'CD5', 'HPGDS', 'IFNG', 'CCR1', 'CD4', 'IL7R', 'LTB'),
  CD8_T = c('CD8A', 'CD8B', 'CD3D', 'CD3E'),
  B_cells = c('MS4A1', 'TCL1A'),
  Treg = c('FOXP3', 'IL2RA', 'TIGIT'),
  Myeloid = c('MMP9', 'LYZ', 'CST3'),
  NK_cells = c('GZMA', 'KLRF1', 'KLRD1', 'CD3D', 'CD3E'),
  Heat_shock = c('HSPA1A', 'HSPA1B', 'DNAJB1')
)

# Generate plots for each marker group
for (group in names(marker_groups)) {
  create_feature_dotplots(
    PBMC_filtered, 
    features = marker_groups[[group]], 
    group.by = 'seurat_clusters', 
    filename = paste0(group, '_markers.pdf')
  )
}

# Batch correction of sample groups with Harmony
PBMC_filtered_harmony <- PBMC_filtered %>% 
  RunHarmony(group.by.vars = c('orig.ident'), plot_convergence = TRUE)

pdf(file = file.path(ann_dir, 'Harmony_batch_correction_plots.pdf'))
DimPlot(PBMC_filtered_harmony, reduction = 'harmony', group.by = 'orig.ident')
VlnPlot(PBMC_filtered_harmony, features = 'harmony_2', group.by = 'orig.ident')
dev.off()

# Clustering the result of Harmony
PBMC_filtered_harmony <- PBMC_filtered_harmony %>%
  RunHarmony(group.by.vars = c('orig.ident'), plot_convergence = TRUE) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(.)) %>%
  RunUMAP(reduction = 'harmony', dims = 1:50, n.neighbors = 15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:50) %>%
  RunTSNE(reduction = 'harmony', dims = 1:50, n.neighbors = 15) %>%
  FindClusters()

# Save Harmony-corrected clustering plots
pdf('PBMC_clustering_after_harmony_seurat.pdf')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'orig.ident')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
DimPlot(PBMC_filtered_harmony, reduction = 'tsne', group.by = 'orig.ident')
DimPlot(PBMC_filtered_harmony, reduction = 'tsne', group.by = 'seurat_clusters')
dev.off()

# Generate plots for each marker group
for (group in names(marker_groups)) {
  create_feature_dotplots(
    PBMC_filtered_harmony, 
    features = marker_groups[[group]], 
    group.by = 'seurat_clusters', 
    filename = paste0(group, '_markers.pdf')
  )
}

# Saving the Harmony corrected output
saveRDS(PBMC_filtered_harmony, 'PBMC_filtered_harmony.rds')

# Save in h5ad format
sceasy::convertFormat(
  PBMC_filtered, 
  from = "seurat", to = "anndata", 
  outFile = 'PBMC_filtered_harmony.h5ad'
)

###################################### Session Info ######################################
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.2

# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# time zone: America/Monterrey
# tzcode source: internal

# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices datasets  utils     methods   base     

# other attached packages:
#  [1] org.Hs.eg.db_3.20.0         AnnotationDbi_1.68.0        clusterProfiler_4.14.4      lubridate_1.9.4             forcats_1.0.0               stringr_1.5.1              
#  [7] purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                tidyverse_2.0.0             SeuratDisk_0.0.0.9021      
# [13] harmony_1.2.3               Rcpp_1.0.14                 cowplot_1.1.3               anndata_0.7.5.6             sceasy_0.0.7                reticulate_1.40.0          
# [19] DoubletFinder_2.0.4         lattice_0.22-6              glmGamPoi_1.18.0            SoupX_1.6.2                 dplyr_1.1.4                 patchwork_1.3.0            
# [25] scran_1.34.0                scater_1.34.0               scuttle_1.16.0              ggplot2_3.5.1               SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
# [31] Biobase_2.66.0              GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0        
# [37] MatrixGenerics_1.18.1       matrixStats_1.5.0           Seurat_4.3.0                SeuratObject_5.0.2          sp_2.2-0                    dior_0.1.5                 

# loaded via a namespace (and not attached):
#   [1] fs_1.6.5                spatstat.sparse_3.1-0   enrichplot_1.26.6       httr_1.4.7              RColorBrewer_1.1-3      tools_4.4.2             sctransform_0.4.1      
#   [8] backports_1.5.0         R6_2.6.0                lazyeval_0.2.2          uwot_0.2.2              withr_3.0.2             gridExtra_2.3           progressr_0.15.1       
#  [15] cli_3.6.4               spatstat.explore_3.3-4  fastDummies_1.7.5       spatstat.data_3.1-4     ggridges_0.5.6          pbapply_1.7-2           yulab.utils_0.2.0      
#  [22] gson_0.1.0              foreign_0.8-88          R.utils_2.12.3          DOSE_4.0.0              parallelly_1.42.0       maps_3.4.2.1            limma_3.62.2           
#  [29] rstudioapi_0.17.1       RSQLite_2.3.9           gridGraphics_0.5-1      generics_0.1.3          ica_1.0-3               spatstat.random_3.3-2   GO.db_3.20.0           
#  [36] Matrix_1.7-2            ggbeeswarm_0.7.2        abind_1.4-8             R.methodsS3_1.8.2       lifecycle_1.0.4         edgeR_4.4.2             qvalue_2.38.0          
#  [43] SparseArray_1.6.1       Rtsne_0.17              blob_1.2.4              grid_4.4.2              promises_1.3.2          dqrng_0.4.1             crayon_1.5.3           
#  [50] ggtangle_0.0.6          miniUI_0.1.1.1          beachmat_2.22.0         KEGGREST_1.46.0         pillar_1.10.1           knitr_1.49              metapod_1.14.0         
#  [57] fgsea_1.32.2            future.apply_1.11.3     codetools_0.2-20        fastmatch_1.1-6         glue_1.8.0              ggfun_0.1.8             spatstat.univar_3.1-1  
#  [64] data.table_1.16.4       treeio_1.30.0           vctrs_0.6.5             png_0.1-8               spam_2.11-1             gtable_0.3.6            assertthat_0.2.1       
#  [71] cachem_1.1.0            xfun_0.50               S4Arrays_1.6.0          mime_0.12               survival_3.8-3          fields_16.3             statmod_1.5.0          
#  [78] bluster_1.16.0          fitdistrplus_1.2-2      ROCR_1.0-11             nlme_3.1-167            ggtree_3.14.0           bit64_4.6.0-1           RcppAnnoy_0.0.22       
#  [85] irlba_2.3.5.1           vipor_0.4.7             KernSmooth_2.23-26      rpart_4.1.24            colorspace_2.1-1        DBI_1.2.3               Hmisc_5.2-2            
#  [92] nnet_7.3-20             tidyselect_1.2.1        bit_4.5.0.1             compiler_4.4.2          htmlTable_2.4.3         BiocNeighbors_2.0.1     hdf5r_1.3.12           
#  [99] DelayedArray_0.32.0     plotly_4.10.4           checkmate_2.3.2         scales_1.3.0            lmtest_0.9-40           digest_0.6.37           goftest_1.2-3          
# [106] spatstat.utils_3.1-2    rmarkdown_2.29          XVector_0.46.0          htmltools_0.5.8.1       pkgconfig_2.0.3         base64enc_0.1-3         fastmap_1.2.0          
# [113] rlang_1.1.5             htmlwidgets_1.6.4       UCSC.utils_1.2.0        shiny_1.10.0            farver_2.1.2            zoo_1.8-12              jsonlite_1.8.9         
# [120] BiocParallel_1.40.0     R.oo_1.27.0             GOSemSim_2.32.0         BiocSingular_1.22.0     magrittr_2.0.3          ggplotify_0.1.2         Formula_1.2-5          
# [127] GenomeInfoDbData_1.2.13 dotCall64_1.2           munsell_0.5.1           ape_5.8-1               viridis_0.6.5           stringi_1.8.4           zlibbioc_1.52.0        
# [134] MASS_7.3-64             plyr_1.8.9              listenv_0.9.1           ggrepel_0.9.6           deldir_2.0-4            Biostrings_2.74.1       splines_4.4.2          
# [141] tensor_1.5              hms_1.1.3               locfit_1.5-9.11         igraph_2.1.4            spatstat.geom_3.3-5     RcppHNSW_0.6.0          reshape2_1.4.4         
# [148] ScaledMatrix_1.14.0     evaluate_1.0.3          renv_1.1.1              BiocManager_1.30.25     tzdb_0.4.0              httpuv_1.6.15           RANN_2.6.2             
# [155] polyclip_1.10-7         future_1.34.0           scattermore_1.2         rsvd_1.0.5              xtable_1.8-4            tidytree_0.4.6          RSpectra_0.16-2        
# [162] later_1.4.1             viridisLite_0.4.2       aplot_0.2.4             memoise_2.0.1           beeswarm_0.4.0          cluster_2.1.8           timechange_0.3.0       
# [169] globals_0.16.3    
