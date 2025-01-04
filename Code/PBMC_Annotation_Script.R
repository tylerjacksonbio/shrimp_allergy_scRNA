##### Annotation using of Whole Data using Azimuth. Gamma-delta T-cell Annotation also Included, Annotated using the Module Score #####
# Importing required libraries
library(Seurat)
library(sceasy)
library(patchwork)
library(cowplot)
library(clustree)
set.seed(42)

# This is the URL for running Azimuth using the web interface - https://app.azimuth.hubmapconsortium.org/app/human-pbmc
# Load the Harmony object and convert formats
PBMC_filtered_harmony <- readRDS('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Datasets/Final_Annotated_Data/PBMC_filtered_harmony.rds')
PBMC_filtered_harmony <- readRDS('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Test_Script/Final_Data_Upload_Analysis/Merged_adata.rds')

# The .h5Seurat format is required for using the Azimuth annotation program
SaveH5Seurat(PBMC_filtered_harmony, 'PBMC_filtered_harmony.h5Seurat')

# Add Azimuth predictions from the annotation program
add_azimuth_metadata <- function(object, filepath) {
  predictions <- read.delim(filepath, row.names = 1)
  AddMetaData(object = object, metadata = predictions)
}

PBMC_filtered_harmony <- add_azimuth_metadata(PBMC_filtered_harmony, 'azimuth_pred_intermediate.tsv')
PBMC_filtered_harmony <- add_azimuth_metadata(PBMC_filtered_harmony, 'azimuth_pred_broad.tsv')

# Add UMAP projections
projected.umap <- readRDS('azimuth_umap.Rds')
PBMC_filtered_harmony <- PBMC_filtered_harmony[, Cells(projected.umap)]
PBMC_filtered_harmony[['umap.proj']] <- projected.umap

# Define plotting function
save_dim_feature_plots <- function(object, filename, reductions, groups, features = NULL, pt.size = 0.5, order = FALSE) {
  pdf(filename)
  for (reduction in reductions) {
    for (group in groups) {
      DimPlot(object, reduction = reduction, group.by = group, pt.size = pt.size, label = FALSE) + theme_minimal()
    }
  }
  if (!is.null(features)) {
    FeaturePlot(object, features, pt.size = pt.size, order = order) + theme_minimal()
  }
  dev.off()
}

# Plot Azimuth classifications and UMAP projections
save_dim_feature_plots(
  PBMC_filtered_harmony, 'Classification_Azimuth.pdf',
  reductions = c('umap.proj', 'umap'),
  groups = c('orig.ident', 'seurat_clusters', 'predicted.celltype.l1', 'predicted.celltype.l2')
)

# Marker gene visualization by cell type after Azimuth annotation
marker_genes <- list(
  "CD4+_Tcell_markers" = c('CSF2', 'CCR10', 'CD5', 'HPGDS', 'IFNG', 'CCR1', 'CD4', 'IL7R', 'LTB'),
  "CD8+_Tcell_markers" = c('CD8A', 'CD8B', 'CD3D', 'CD3E'),
  "B-cell_markers" = c('MS4A1', 'TCL1A'),
  "Treg_markers" = c('FOXP3', 'IL2RA', 'TIGIT'),
  "Myeloid_markers" = c('MMP9', 'LYZ', 'CST3'),
  "NK_cell_markers" = c('GZMA', 'KLRF1', 'KLRD1', 'CD3D', 'CD3E'),
  "Heat_shock_rich" = c('HSPA1A', 'HSPA1B', 'DNAJB1'),
  "gdT_cell_markers" = c('TRDC', 'TRDV1', 'TRDV2', 'NKG7', 'CD3D', 'CD3E', 'CD8A', 'CD8B', 'CD4')
)

for (name in names(marker_genes)) {
  pdf(paste0(name, '.pdf'))
  FeaturePlot(PBMC_filtered_harmony, marker_genes[[name]], pt.size = 0.5, order = TRUE) + theme_minimal()
  plot_grid(
    DotPlot(
      object = PBMC_filtered_harmony,
      features = marker_genes[[name]],
      group.by = 'seurat_clusters'
    ) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
  )
  dev.off()
}

# Set final annotations and plot all markers
Idents(PBMC_filtered_harmony) <- 'Final_annotation'
pdf('All_markers.pdf', height = 10, width = 10)
plot_grid(
  DotPlot(
    object = PBMC_filtered_harmony,
    features = unlist(marker_genes),
    group.by = 'Final_annotation'
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
)
dev.off()

##### Reannotating the data using Azimuth predictions and marker genes #####
# Find clusters at the desired resolution and set up final annotation metadata
PBMC_filtered_harmony <- readRDS('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Test_Script/Final_Data_Upload_Analysis/Merged_adata.rds')
PBMC_filtered_harmony<-FindNeighbors(PBMC_filtered_harmony, reduction = 'harmony', dims = 1:50)
PBMC_filtered_harmony<-RunUMAP(PBMC_filtered_harmony, reduction = 'harmony', dims = 1:50, n.neighbors = 15)
PBMC_filtered_harmony<-RunTSNE(PBMC_filtered_harmony, reduction = 'harmony', dims = 1:50, n.neighbors = 15)
resolutions <- c(2, 1, 0.8, 0.5, 0.4, 0.3, 0.2, 0.1)

# Apply FindClusters for each resolution
for (res in resolutions) {
  PBMC_filtered_harmony <- FindClusters(PBMC_filtered_harmony, resolution = res)
}

pdf('PBMC_filtered_harmony_clustering.pdf')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'seurat_clusters')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE) + NoLegend()
DimPlot(PBMC_filtered_harmony, reduction = 'tsne', group.by = 'seurat_clusters')
dev.off()

# Set identity and build cluster tree
Idents(PBMC_filtered_harmony) <- 'Final_annotation_broad'
PBMC_filtered_harmony <- BuildClusterTree(PBMC_filtered_harmony)

# Extract phylogenetic tree
myPhyTree <- Tool(object = PBMC_filtered_harmony, slot = "BuildClusterTree")

# Save plots to PDF
pdf('ClusterTrees_Whole_Data.pdf', height = 12, width = 8)
clustree(PBMC_filtered_harmony, prefix = "RNA_snn_res.")
PlotClusterTree(PBMC_filtered_harmony, direction = "rightwards")

pdf('ClusterTrees_Whole_Data_ggtree.pdf')
ggtree(myPhyTree) + geom_tiplab() + theme_tree() + xlim(NA, 1000)
dev.off()

# Select resolution of 0.5 for annotation (best res based on dendrogram)
Idents(PBMC_filtered_harmony) <- 'RNA_snn_res.0.5'
PBMC_filtered_harmony[["Final_annotation_broad"]] <- Idents(PBMC_filtered_harmony)

# Rename cluster identities with clear labels for final annotation
cluster_annotations <- c(
  "0" = "CD4 T-cells", "1" = "CD4 T-cells", "2" = "CD8 T-cells",
  "3" = "CD8 T-cells", "4" = "Naive B-cells", "5" = "NK cells",
  "6" = "Memory/Intermediate B-cells", "7" = "CD8 T-cells", "8" = "Broad cell types",
  "9" = "Broad cell types", "10" = "Treg cells", "11" = "Broad cell types",
  "12" = "Monocytes", "13" = "CD8 T-cells", "14" = "CD4 T-cells",
  "15" = "NK cells"
)
PBMC_filtered_harmony <- RenameIdents(PBMC_filtered_harmony, cluster_annotations)

# Intermediate annotations
PBMC_filtered_harmony[["Final_annotation_intermediate"]] <- Idents(object = PBMC_filtered_harmony)
Idents(PBMC_filtered_harmony) <- PBMC_filtered_harmony@meta.data$predicted.celltype.l2
PBMC_filtered_harmony <- RenameIdents(object = PBMC_filtered_harmony, `B intermediate` = "B intermediate", `B memory` = "B memory", `B naive` = "B naive",
                                      `CD14 Mono` = 'CD14 Mono',`CD16 Mono` = "CD16 Mono", `CD4 CTL` = "CD4 CTL", `CD4 Proliferating` = "CD4 Proliferating",
                                      `CD4 TCM` = 'CD4 TCM',`CD4 TEM` = "CD4 TEM", `CD8 Naive` = "CD8 Naive", `CD8 TCM` = "CD8 TCM",
                                      `CD8 TEM` = 'CD8 TEM',`cDC2` = "cDC2", `dnT` = "dnT", `Eryth` = "Eryth",
                                      `gdT` = "Reann", `dnT` = "dnT", `HSPC` = "HSPC", `ILC` = "ILC", `MAIT` = "MAIT",
                                      `NK` = "NK", `NK Proliferating` = "NK Proliferating", `NK_CD56bright` = "NK_CD56bright", `Plasmablast` = "Plasmablast", `Platelet` = "Platelet",
                                      `Treg` = "Treg cells")

# Plot new annotations
plot_annotations <- function(object, filename) {
  pdf(filename, height = 12, width = 10)
  DimPlot(object, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 1, label = TRUE) +
    ggtitle("UMAP: Original Seurat Clusters") + theme_minimal()
  DimPlot(object, reduction = 'umap', group.by = 'Final_annotation', pt.size = 1) +
    ggtitle("UMAP: Final Annotation (No Labels)") + theme_minimal()
  DimPlot(object, reduction = 'umap', group.by = 'Final_annotation', label = TRUE, label.size = 6, pt.size = 1) +
    ggtitle("UMAP: Final Annotation (With Labels)") + theme_minimal()
  dev.off()
}

plot_annotations(PBMC_filtered_harmony, "Final_annotation.pdf")

# Save the annotated object and convert format
output_path <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Annotated_Data/"
saveRDS(PBMC_filtered_harmony, file.path(output_path, "Final_annotated_object.rds"))
sceasy::convertFormat(
  PBMC_filtered_harmony,
  from = "seurat",
  to = "anndata",
  outFile = file.path(output_path, "Final_annotated_object.h5ad")
)

##### Save the final annotated object used for downstream analysis #####
saveRDS(PBMC_filtered_harmony, file.path(output_path, "Merged_adata.rds"))
sceasy::convertFormat(
  PBMC_filtered_harmony,
  from = "seurat",
  to = "anndata",
  outFile = file.path(output_path, "Merged_adata.h5ad")
)

##### Whole Data Cell Proportion Analysis and Preliminary Gamma-delta T-cell Analysis #####
##### Whole data cell proportion analysis #####
# Set working directory and define color palette
setwd('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/gdT_cell_annotation')
sample_colors <- c('HC-unstim' = '#9FBDD3', 'HC-stim' = '#5A6EDC', 
                   'SA-unstim' = '#F6CC93', 'SA-stim' = '#F09924')

# Assign identities
Idents(PBMC_filtered_harmony) <- PBMC_filtered_harmony@meta.data$orig.ident

# Subset samples and create tables
samples <- list(
  HC_unstim = subset(PBMC_filtered_harmony, orig.ident == 'HC-unstim'),
  HC_stim = subset(PBMC_filtered_harmony, orig.ident == 'HC-stim'),
  SA_unstim = subset(PBMC_filtered_harmony, orig.ident == 'SA-unstim'),
  SA_stim = subset(PBMC_filtered_harmony, orig.ident == 'SA-stim')
)

# Output annotation tables for each sample
lapply(samples, function(sample) table(sample@meta.data$Final_annotation_broad))

# Plot UMAPs for all samples
pdf('UMAPS_split_by_sample.pdf')
for (sample_name in names(samples)) {
  DimPlot(samples[[sample_name]], reduction = 'umap', group.by = 'orig.ident', 
          cols = sample_colors[sample_name], pt.size = 0.1) +
    NoLegend() + ggtitle(sample_name)
}
dev.off()

##### Annotating and reclustering the gdT-cells #####
# Annotating gamma-delta T-cells
PBMC_filtered_harmony$orig.ident <- factor(PBMC_filtered_harmony$orig.ident, 
                                           levels = names(sample_colors))
PBMC_filtered_harmony$group.colors <- sample_colors[as.character(PBMC_filtered_harmony$orig.ident)]

gdt_markers <- c('TRDC','CD7','CD247','GZMA','SPON2','CTSW','KLRD1','CST7','HOPX','PRF1','FGFBP2','GZMB','CCL5','CD3D','CD3E')
PBMC_filtered_harmony <- AddModuleScore(PBMC_filtered_harmony, features = gdt_markers, ctrl = 5, name = 'gdt_Module_Score')
PBMC_filtered_harmony <- FindClusters(PBMC_filtered_harmony, resolution = 2)

# Plot the result of the module score
pdf('module_score_gdt.pdf')
DimPlot(PBMC_filtered_harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
FeaturePlot(PBMC_filtered_harmony, features = 'gdt_Module_Score1', cols = c('lightgrey', 'red'), order = TRUE , pt.size = .7)
VlnPlot(PBMC_filtered_harmony, features = 'gdt_Module_Score1', assay = 'RNA')
dev.off()

# Subset the cells that show a cutoff of 0.5 for the module score
Idents(PBMC_filtered_harmony) <- PBMC_filtered_harmony@meta.data$seurat_clusters
Gdt_cells <- subset(PBMC_filtered_harmony, subset = seurat_clusters ==  '11' | seurat_clusters == '17' 
                    | seurat_clusters == '20')

pdf('Gdt_cell_cluster_subset_before_recluster.pdf')
DimPlot(Gdt_cells, reduction = 'umap', group.by = 'seurat_clusters')
VlnPlot(Gdt_cells, features = 'gdt_Module_Score1', assay = 'RNA')
dev.off()

Gdt_cells

Gdt_cells <- subset(Gdt_cells, subset = gdt_Module_Score1 > 0.5)
table(Gdt_cells@meta.data$orig.ident)

# We will annotate the clusters that are above the module score cutoff of 0.5 as gdt cells
# Set annotation column in gdt_cells to "Gamma_delta_T-cells"
Gdt_cells@meta.data$Final_annotation <- "gdT-cells"

# Create a named vector of annotation column from gdt_cells
dict_D <- setNames(Gdt_cell@meta.data$Final_annotation_broad, colnames(Gdt_cell_clusters))

# Get a list of Final_annotation from PBMC_Data, check the barcode in dict_D
new_L <- character(length(PBMC_filtered_harmony))
new_L
for (i in 1:ncol(PBMC_filtered_harmony)){
  barcode <- colnames(PBMC_filtered_harmony)[i]
  annot <- PBMC_filtered_harmony@meta.data$Final_annotation_broad[i]
  if (barcode %in% names(dict_D)) new_L[i] <- dict_D[barcode] else new_L[i] <- annot
}

# Set New_Annotation column in PBMC_Data to new_L
PBMC_filtered_harmony@meta.data$Final_annotation_broad <- factor(new_L)
table(PBMC_filtered_harmony@meta.data$Final_annotation_broad)

# Subset gamma-delta T-cells
Gdt_cell_clusters <- subset(PBMC_filtered_harmony, Final_annotation == 'gd T-cells')

# Visualize gamma-delta T-cells
pdf('gdT_cells_subset.pdf')
DimPlot(PBMC_filtered_harmony, group.by = 'Final_annotation')
dev.off()

pdf('Gdt_cells_before_reclustering.pdf')
DimPlot(Gdt_cell_clusters, reduction = 'umap', group.by = 'seurat_clusters') +
  DimPlot(Gdt_cell_clusters, reduction = 'umap', group.by = 'orig.ident', cols = sample_colors)
dev.off()

# Checking the appropriate number of dimensions
PBMC_filtered_harmony <- readRDS('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Datasets/Final_Annotated_Data/Final_Data_Upload_Analysis/Merged_adata.rds')

# Set the identity class and subset data
Idents(PBMC_filtered_harmony) <- 'Final_annotation_broad'
Gdt_cell_clusters <- subset(PBMC_filtered_harmony, subset = Final_annotation_broad == 'gd T-cells')

# Define a function to process clustering
process_clustering <- function(object, dims, resolution = 0.3, neighbors = 15) {
  object <- FindNeighbors(object, reduction = 'harmony', dims = 1:dims)
  ElbowPlot(object, ndims = dims)
  object <- RunUMAP(object, reduction = 'harmony', dims = 1:dims, n.neighbors = neighbors)
  object <- RunTSNE(object, reduction = 'harmony', dims = 1:dims, n.neighbors = neighbors)
  object <- FindClusters(object, resolution = resolution)
  return(object)
}

# Loop through dimensions and process clustering
dims_list <- c(50, 40, 30, 20, 10)
for (dims in dims_list) {
  Gdt_cell_clusters <- process_clustering(Gdt_cell_clusters, dims = dims)
  
  # Generate plots and save to PDF
  pdf(paste0(dims, '_dims_gdt_cell_clustering.pdf'))
  DimPlot(Gdt_cell_clusters, group.by = 'seurat_clusters', reduction = 'umap') + 
    ggtitle('Seurat Clusters')
  DimPlot(Gdt_cell_clusters, group.by = 'Gdt_cluster_subsets', reduction = 'umap') + 
    ggtitle('gdT Subclusters')
  DimPlot(Gdt_cell_clusters, group.by = 'orig.ident', cols = sample_colors, reduction = 'umap') + 
    ggtitle('Sample Group')
  DimPlot(Gdt_cell_clusters, group.by = 'seurat_clusters', reduction = 'umap', label = FALSE) + 
    ggtitle('Seurat Clusters') + NoLegend()
  DimPlot(Gdt_cell_clusters, group.by = 'Gdt_cluster_subsets', reduction = 'umap', label = FALSE) + 
    ggtitle('gdT Subclusters') + NoLegend()
  DimPlot(Gdt_cell_clusters, group.by = 'orig.ident', cols = sample_colors, reduction = 'umap', label = FALSE) + 
    ggtitle('Sample Group') + NoLegend()
  dev.off()
}

# Reannotate the seurat clusters
Gdt_cell_clusters[["Gdt_cluster_subsets"]] <- Idents(Gdt_cell_clusters)
Idents(Gdt_cell_clusters) <- Gdt_cell_clusters$seurat_clusters
Gdt_cell_clusters <- RenameIdents(Gdt_cell_clusters, 
                                  `0` = "gdt_cluster_1", 
                                  `1` = "gdt_cluster_2", 
                                  `2` = "gdt_cluster_3")

# Building a dendrogram for clustering in the gdT-cells
Gdt_cell_clusters <- FindNeighbors(Gdt_cell_clusters, reduction = 'harmony', dims = 1:50)
ElbowPlot(Gdt_cell_clusters, ndims = 50) 
Gdt_cell_clusters <- RunUMAP(Gdt_cell_clusters, reduction = 'harmony', dims = 1:50, n.neighbors = 15) 
Gdt_cell_clusters <- RunTSNE(Gdt_cell_clusters, reduction = 'harmony', dims = 1:50, n.neighbors = 15) 

# Define resolutions to test
resolutions <- c(2, 1, 0.8, 0.5, 0.4, 0.3, 0.2, 0.1)

# Apply FindClusters for each resolution
for (res in resolutions) {
  Gdt_cell_clusters <- FindClusters(Gdt_cell_clusters, resolution = res)
}

# Set identity and build cluster tree
Idents(Gdt_cell_clusters) <- 'Gdt_cluster_subsets'
Gdt_cell_clusters <- BuildClusterTree(Gdt_cell_clusters)

# Extract phylogenetic tree
myPhyTree <- Tool(object = Gdt_cell_clusters, slot = "BuildClusterTree")

# Save plots to PDF
pdf('ClusterTrees_gdT_cells.pdf', height = 12, width = 8)
clustree(Gdt_cell_clusters, prefix = "RNA_snn_res.")
PlotClusterTree(Gdt_cell_clusters, direction = "rightwards")
ggtree(myPhyTree) + geom_tiplab() + theme_tree() + xlim(NA, 10)
dev.off()

# Rename and annotate clusters
Gdt_cell_clusters[["Gdt_cluster_subsets"]] <- Idents(Gdt_cell_clusters)
Idents(Gdt_cell_clusters) <- Gdt_cell_clusters$seurat_clusters
Gdt_cell_clusters <- RenameIdents(Gdt_cell_clusters, `0` = "gdt_cluster_1", 
                                  `1` = "gdt_cluster_2", `2` = "gdt_cluster_3")

# Visualize gamma-delta T-cell clusters
pdf('Module_scores_and_gdt_cell_subsets.pdf')
DimPlot(Gdt_cell_clusters, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE) +
  DimPlot(Gdt_cell_clusters, reduction = 'umap', group.by = 'Gdt_cluster_subsets', label = TRUE) +
  FeaturePlot(Gdt_cell_clusters, features = 'gdt_Module_Score1', cols = c('lightgrey', 'red'), 
              order = TRUE, pt.size = 0.7) +
  VlnPlot(Gdt_cell_clusters, features = 'gdt_Module_Score1', assay = 'RNA')
dev.off()

