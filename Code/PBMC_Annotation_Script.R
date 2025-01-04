##### Annotation using of Whole Data using Azimuth. Gamm-delta T-cell Annotation also Included, Annotated using the Module Score #####
# Importing required libraries
library(Seurat)
library(sceasy)
library(patchwork)
library(cowplot)

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
PBMC_filtered_harmony <- FindClusters(PBMC_filtered_harmony, resolution = 0.5)
PBMC_filtered_harmony[["Final_annotation"]] <- Idents(PBMC_filtered_harmony)

# Rename cluster identities with clear labels for final annotation
cluster_annotations <- c(
  "0" = "CD4_T_cells", "1" = "CD4_T_cells", "2" = "CD8_T_cells",
  "3" = "CD8_T_cells", "4" = "Naive_B_cells", "5" = "NK_cells",
  "6" = "Memory_B_cells", "7" = "CD8_T_cells", "8" = "Other_T_cells",
  "9" = "Other_T_cells", "10" = "Treg_cells", "11" = "Other_T_cells",
  "12" = "Monocytes", "13" = "CD8_T_cells", "14" = "CD4_T_cells",
  "15" = "NK_cells"
)
PBMC_filtered_harmony <- RenameIdents(PBMC_filtered_harmony, cluster_annotations)

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
output_path <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/With_ambient_removed/"
saveRDS(PBMC_filtered_harmony, file.path(output_path, "Final_annotated_object.rds"))
sceasy::convertFormat(
  PBMC_filtered_harmony,
  from = "seurat",
  to = "anndata",
  outFile = file.path(output_path, "Final_annotated_object.h5ad")
)

# Define constants
output_path <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/With_ambient_removed/"
gdt_markers <- c('TRDC', 'CD7', 'CD247', 'GZMA', 'SPON2', 'CTSW', 'KLRD1', 
                 'CST7', 'HOPX', 'PRF1', 'FGFBP2', 'GZMB', 'CCL5', 'CD3D', 'CD3E')

# Add module score for gamma-delta T-cells and recluster
PBMC_filtered_harmony <- AddModuleScore(PBMC_filtered_harmony, features = gdt_markers, ctrl = 5, name = "gdt_Module_Score")
PBMC_filtered_harmony <- FindClusters(PBMC_filtered_harmony, resolution = 2)

# Plot module score results
pdf(file.path(output_path, "module_score_gdt.pdf"))
DimPlot(PBMC_filtered_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
FeaturePlot(PBMC_filtered_harmony, features = "gdt_Module_Score1", cols = c("lightgrey", "red"), order = TRUE, pt.size = 0.7)
VlnPlot(PBMC_filtered_harmony, features = "gdt_Module_Score1", assay = "RNA")
dev.off()

# Subset gamma-delta T-cells based on clusters and module score cutoff
gdt_clusters <- c("11", "17", "20")
Gdt_cells <- subset(PBMC_filtered_harmony, seurat_clusters %in% gdt_clusters & gdt_Module_Score1 > 0.5)

# Plot subset results before reclustering
pdf(file.path(output_path, "Gdt_cell_cluster_subset_before_recluster.pdf"))
DimPlot(Gdt_cells, reduction = "umap", group.by = "seurat_clusters")
VlnPlot(Gdt_cells, features = "gdt_Module_Score1", assay = "RNA")
dev.off()

# Re-run dimensionality reduction and clustering
Gdt_cells <- RunUMAP(Gdt_cells, reduction = "harmony", dims = 1:50, n.neighbors = 15)
Gdt_cells <- FindNeighbors(Gdt_cells, reduction = "harmony", dims = 1:50)
Gdt_cells <- FindClusters(Gdt_cells, resolution = 0.3)

# Plot subset results after reclustering
pdf(file.path(output_path, "Gdt_cell_cluster_subset_after_recluster.pdf"))
DimPlot(Gdt_cells, reduction = "umap", group.by = "seurat_clusters")
dev.off()

# Annotate gamma-delta T-cells in the main dataset
PBMC_filtered_harmony$Final_annotation <- as.character(PBMC_filtered_harmony$Final_annotation)
gamma_delta_annotation <- rep("Gamma_delta_T-cells", length(colnames(Gdt_cells)))
names(gamma_delta_annotation) <- colnames(Gdt_cells)
PBMC_filtered_harmony$Final_annotation[colnames(Gdt_cells)] <- gamma_delta_annotation
PBMC_filtered_harmony$Final_annotation <- factor(PBMC_filtered_harmony$Final_annotation)

# Final plot for annotation
pdf(file.path(output_path, "Final_annotation_NEW.pdf"))
DimPlot(PBMC_filtered_harmony, group.by = "Final_annotation", reduction = "umap", label = TRUE)
DimPlot(PBMC_filtered_harmony, group.by = "Final_annotation", reduction = "tsne")
dev.off()

##### Save the final annotated object used for downstream analysis #####
saveRDS(PBMC_filtered_harmony, file.path(output_path, "Merged_adata.rds"))
sceasy::convertFormat(
  PBMC_filtered_harmony,
  from = "seurat",
  to = "anndata",
  outFile = file.path(output_path, "Merged_adata.h5ad")
)
