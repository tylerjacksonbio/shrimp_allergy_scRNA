##### Importing and Preprocessing Script
##### Loading in the required packages #####
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

##### Performing SoupX Ambient RNA Removal #####
# Loading in the .h5 outputs from CellRanger 6.1.2
setwd('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Cellranger_Output_Files/')
SD1_counts <- Read10X_h5("SD1_filtered_feature_bc_matrix.h5")
SD1_raw_matrix  <- Read10X_h5("SD1_raw_feature_bc_matrix.h5",use.names = T)
SD2_counts <- Read10X_h5("SD2_filtered_feature_bc_matrix.h5")
SD2_raw_matrix  <- Read10X_h5("SD2_raw_feature_bc_matrix.h5",use.names = T)
SD3_counts <- Read10X_h5("SD3_filtered_feature_bc_matrix.h5")
SD3_raw_matrix  <- Read10X_h5("SD3_raw_feature_bc_matrix.h5",use.names = T)
SD4_counts <- Read10X_h5("SD4_filtered_feature_bc_matrix.h5")
SD4_raw_matrix  <- Read10X_h5("SD4_raw_feature_bc_matrix.h5",use.names = T)

# Define datasets and names
datasets <- list(SD1 = list(counts = SD1_counts, raw_matrix = SD1_raw_matrix),
                 SD2 = list(counts = SD2_counts, raw_matrix = SD2_raw_matrix),
                 SD3 = list(counts = SD3_counts, raw_matrix = SD3_raw_matrix),
                 SD4 = list(counts = SD4_counts, raw_matrix = SD4_raw_matrix))
dataset_names <- names(datasets)

# Base output directory
output_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/Correct_PBMC_Analysis/QC_Figures/SoupX_Figures/Before_Removal/"

output_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Test_Script/Importing_Preprocessing"

srat <- CreateSeuratObject(counts = datasets$SD1$counts)
soup.channel <- SoupChannel(srat$raw_matrix, srat$counts)

# Loop over each dataset to view clustering prior to SoupX
for (i in seq_along(datasets)) {
  dataset_name <- dataset_names[i]
  dataset <- datasets[[dataset_name]]
  
  # Create Seurat object
  srat <- CreateSeuratObject(counts = dataset$counts)
  
  # Create SoupX object
  soup.channel <- SoupChannel(dataset$raw_matrix, dataset$counts)
  
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
  
  # Set output directory for this dataset
  setwd(output_dir)
  
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

# Define datasets and names
datasets <- list(SD1 = list(counts = SD1_counts, raw_matrix = SD1_raw_matrix, out_folder = "soupXOutFolder_SD1/"),
                 SD2 = list(counts = SD2_counts, raw_matrix = SD2_raw_matrix, out_folder = "soupXOutFolder_SD2/"),
                 SD3 = list(counts = SD3_counts, raw_matrix = SD3_raw_matrix, out_folder = "soupXOutFolder_SD3/"),
                 SD4 = list(counts = SD4_counts, raw_matrix = SD4_raw_matrix, out_folder = "soupXOutFolder_SD4/"))
dataset_names <- names(datasets)

# Loop through datasets to run the SoupX algorithm
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
  DropletUtils:::write10xCounts(dataset$out_folder, adj.matrix)
  sratSoupx <- CreateSeuratObject(counts = adj.matrix)
  
  # QC metrics comparison
  srat@meta.data[] <- "Before_SoupX"
  metadata <- srat@meta.data
  colnames(metadata)[] <- "Ambient_RNA_Removal"
  srat@meta.data <- metadata
  
  sratSoupx@meta.data[1] <- "After_SoupX"
  metadata <- sratSoupx@meta.data
  colnames(metadata)[1] <- "Ambient_RNA_Removal"
  sratSoupx@meta.data <- metadata
  
  Merge <- merge(x = srat, y = sratSoupx, add.cell.ids = c("noAmb", "Amb"), project = dataset_name)
  Merge[['percent.mt']] <- PercentageFeatureSet(Merge, pattern = "^MT-")
  Idents(Merge) <- Merge@meta.data$Ambient_RNA_Removal
  Merge@meta.data$Ambient_RNA_Removal <- factor(Merge@meta.data$Ambient_RNA_Removal, levels = c("Before_SoupX", "After_SoupX"))
  
  # Save QC plots
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

# Save outputs for each sample group
for (i in seq_along(datasets)) {
  dataset_name <- dataset_names[i]
  output_folder <- datasets[[dataset_name]]$out_folder
  
  # Read saved data
  adj_data <- Read10X(output_folder)
  srat <- CreateSeuratObject(counts = adj_data, project = dataset_name)
  
  saveRDS(srat, file = paste0(dataset_name, "_srat.rds"))
}

##### Performing Doublet Removal using DoubletFinder #####

# Define datasets and their orig.ident values
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

pdf('QC_plots_after_SoupX.pdf')
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
saveRDS(all.samples, 'Combined_Shrimp_Allergy_before_filtering_and_DoubletFinder.rds')

# Read back the saved file if needed
all.samples <- readRDS('Combined_Shrimp_Allergy_before_filtering_and_DoubletFinder.rds')

# Run the DoubletFinder algorithm to remove doublets from all four samples using for loop
all.samples.split <- SplitObject(all.samples, split.by = "orig.ident") 
for (i in 1:length(all.samples.split)) {
  print(paste0("orig.ident",i))
  PBMC_sample    <- NormalizeData(all.samples.split[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  PBMC_sample    <- FindVariableFeatures(PBMC_sample, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(PBMC_sample)
  PBMC_sample    <- ScaleData(PBMC_sample, features = all.genes)
  PBMC_sample    <- RunPCA(PBMC_sample, npcs = 50, verbose = F)
  #
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
  #
  PBMC_sample <- RunUMAP(PBMC_sample, dims = 1:min.pc)
  PBMC_sample <- FindNeighbors(object = PBMC_sample, dims = 1:min.pc)              
  PBMC_sample <- FindClusters(object = PBMC_sample)
  #
  sweep.list <- paramSweep_v3(PBMC_sample, PCs = 1:min.pc, sct = F)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  #
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  #
  annotations <- PBMC_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(0.035*nrow(PBMC_sample@meta.data)) ## Assuming 3.5% according to 10X recommendations
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  #
  #
  PBMC_sample <- doubletFinder_v3(seu = PBMC_sample, 
                                  PCs = 1:min.pc, 
                                  pK = optimal.pk,
                                  nExp = nExp.poi.adj,
                                  sct = F)
  metadata <- PBMC_sample@meta.data
  colnames(metadata)[8] <- "doublet"
  PBMC_sample@meta.data <- metadata 
  
  # subset and save
  PBMC_sample.singlets <- subset(PBMC_sample, doublet == "Singlet")
  all.samples.split[[i]] <- PBMC_sample.singlets
  remove(PBMC_sample.singlets)
}

table(PBMC_sample@meta.data[["doublet"]])

##### Merge Samples into a New List Containing Only Singlets #####
PBMC_sample.singlets <- merge(
  x = all.samples.split[[1]],
  y = all.samples.split[-1],
  project = "shrimp_pbmc"
)
print(PBMC_sample.singlets)

# Save and reload the singlets data
saveRDS(PBMC_sample.singlets, 'Combined_after_ambient_and_doublet_removal.rds')
PBMC_sample.singlets <- readRDS('Combined_after_ambient_and_doublet_removal.rds')

##### Cluster Data to Compare Before and After Doublet Removal #####
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

##### Quality Control and Filtering #####
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

##### Filter low-quality cells #####
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

##### Apply Harmony Batch Correction for the Different Sample Groups #####
# Plotting the results of the data before batch correction
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

# Directory to save the plots
setwd('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/Correct_PBMC_Analysis/Annotation_Plots/')

pdf('Clustering_after_QC.pdf')
DimPlot(PBMC_filtered, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 0.1, label = TRUE)
DimPlot(PBMC_filtered, reduction = 'umap', group.by = 'orig.ident', cols = my_colors, pt.size = 0.1)
DimPlot(PBMC_filtered, reduction = 'tsne', group.by = 'seurat_clusters', pt.size = 0.1)
dev.off()

# Define a helper function
create_feature_dotplots <- function(data, features, group.by, filename, point.size = 0.7) {
  pdf(filename)
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

# Applying batch correction with Harmony
PBMC_filtered_harmony <- PBMC_filtered %>%
  RunHarmony(group.by.vars = c('orig.ident'), plot_convergence = TRUE) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(.)) %>%
  RunUMAP(reduction = 'harmony', dims = 1:50, n.neighbors = 15) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:50) %>%
  RunTSNE(reduction = 'harmony', dims = 1:50, n.neighbors = 15) %>%
  FindClusters()

# Save Harmony-corrected clustering results
pdf('PBMC_clustering_after_harmony_seurat.pdf')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'orig.ident')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
DimPlot(PBMC_filtered_harmony, reduction = 'tsne', group.by = 'orig.ident')
DimPlot(PBMC_filtered_harmony, reduction = 'tsne', group.by = 'seurat_clusters')
dev.off()

# Saving the Harmony corrected output
saveRDS(PBMC_filtered_harmony, 'PBMC_filtered_harmony.rds')

# Save in h5ad format
sceasy::convertFormat(
  PBMC_filtered, 
  from = "seurat", to = "anndata", 
  outFile = 'PBMC_filtered_harmony.h5ad'
)
