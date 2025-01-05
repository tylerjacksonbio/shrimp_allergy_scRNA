##### This is the Script for Running CellChat #####
# Load libraries
suppressMessages({
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(future)
  library(org.Hs.eg.db)
  library(reticulate)
})

# Set global options
options(stringsAsFactors = FALSE)
set.seed(42)

# Defining the needed directories and paths
base_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su"
data_dir <- file.path(base_dir, "PBMC_Dataset/Final_Datasets/Final_Annotated_Data")
plot_dir <- file.path(base_dir, "PBMC_Plots/Correct_PBMC_Analysis/Cellchat_Analysis/Comparison_Only_LR_Pairs")
save_dir <- file.path(data_dir, "Cellchat_Analysis_Datasets_Only_LR_Pairs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Load dataset
PBMC_filtered_harmony <- readRDS(file.path(data_dir, "Final_Annotated_Data.rds"))
Idents(PBMC_filtered_harmony) <- "Final_annotation_broad"

# Preparing data for CellChat
prepare_cellchat_data <- function(object, condition) {
  subset_obj <- subset(object, subset = orig.ident == condition)
  data_input <- GetAssayData(subset_obj, assay = "RNA", slot = "data")
  meta <- data.frame(group = Idents(subset_obj), row.names = names(Idents(subset_obj)))
  return(list(data = data_input, meta = meta))
}

# The conditions that will be compared
conditions <- c("HC-unstim", "HC-stim", "SA-unstim", "SA-stim")

# Creating the CellChat objects
create_and_save_cellchat <- function(condition) {
  prepared_data <- prepare_cellchat_data(PBMC_filtered_harmony, condition)
  cellchat <- createCellChat(object = prepared_data$data, meta = prepared_data$meta, group.by = "group")
  saveRDS(cellchat, file = file.path(save_dir, paste0("cellchat_", condition, "_data.rds")))
  return(cellchat)
}

cellchat_objects <- map(conditions, create_and_save_cellchat)

# Make function for running CellChat on individual conditions
process_cellchat <- function(cellchat_obj) {
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  cellchat_obj <- computeCommunProb(cellchat_obj)
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  cellchat_obj <- aggregateNet(cellchat_obj)
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  return(cellchat_obj)
}

# Run CellChat processing pipeline, can set the multithreading options if parallelization desired
future::plan("multisession", workers = 4)
processed_cellchat_objects <- map(cellchat_objects, process_cellchat)

# Save processed CellChat objects
walk2(processed_cellchat_objects, conditions, ~ saveRDS(.x, file.path(save_dir, paste0("cellchat_", .y, "_processed.rds"))))

# Comparison analysis
compare_conditions <- list(
  HC_vs_SA = list(HC = processed_cellchat_objects[[1]], SA = processed_cellchat_objects[[3]]),
  HCstim_vs_SAstim = list(HC_stim = processed_cellchat_objects[[2]], SA_stim = processed_cellchat_objects[[4]])
)

perform_comparison <- function(comparison_list, comparison_name) {
  cellchat <- mergeCellChat(comparison_list, add.names = names(comparison_list))
  saveRDS(cellchat, file = file.path(save_dir, paste0("cellchat_comparison_", comparison_name, ".rds")))
}

walk2(compare_conditions, names(compare_conditions), perform_comparison)


# Set the working directory for each comparison
comparisons <- c("HCstim_vs_HCunstim", "SAstim_vs_SAunstim", 
                 "SAstim_vs_HCstim", "SAunstim_vs_HCunstim")

base_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/Correct_PBMC_Analysis/Cellchat_Analysis/Comparison_Only_LR_Pairs"

# Create directories for final plots
lapply(comparisons, function(comp) {
  dir.create(file.path(base_dir, comp, "Final_Plots"), showWarnings = FALSE, recursive = TRUE)
})

# Iterate through each comparison
for (comp in comparisons) {
  comp_dir <- file.path(base_dir, comp)
  setwd(comp_dir)
  
  # Compare the number of interactions and the interaction strength
  gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
  gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
  
  pdf(file.path(comp_dir, "Differential_interactions.pdf"))
  print(gg1 + gg2)
  dev.off()
  
  # Visualize differential interactions with heatmaps
  heatmap1 <- netVisual_heatmap(cellchat)
  heatmap2 <- netVisual_heatmap(cellchat, measure = "weight")
  
  pdf(file.path(comp_dir, "Heatmap_interactions.pdf"))
  print(heatmap1 + heatmap2)
  dev.off()
  
  # Visualize interactions between cell populations
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, arrow.width = 2, arrow.size = 0.5, edge.width.max = 5)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight", arrow.width = 2, edge.width.max = 8, targets.use = "gd T-cells")
  
  # Set dataset factors and plot differentially expressed genes
  cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = strsplit(comp, "_vs_")[[1]])
  signaling_pathways <- list(
    "HCstim_vs_HCunstim" = c("MHC-II", "ITGB2", "MHC-I"),
    "SAstim_vs_SAunstim" = c("MHC-II", "ITGB2", "CLEC", "SELPLG", "SPP1"),
    "SAunstim_vs_HCunstim" = c("PARs", "TGFb", "ITGB2", "MHC-II"),
    "SAstim_vs_HCstim" = c("TGFb", "BTLA", "TNF", "CLEC")
  )
  
  pdf(file.path(comp_dir, "Signaling_pathway_gene_expression.pdf"))
  for (pathway in signaling_pathways[[comp]]) {
    plotGeneExpression(cellchat, signaling = pathway, split.by = "datasets", colors.ggplot = TRUE)
  }
  dev.off()
  
  # Save CellChat object for the comparison
  saveRDS(cellchat, file.path(base_dir, "Final_Annotated_Data", paste0("cellchat_", comp, "_comparison.rds")))
}

# Plot final signaling pathways in gdT cells
gdT_cell_signaling <- list(
  "HCstim_vs_HCunstim" = c("IL16", "TGFb"),
  "SAstim_vs_SAunstim" = c("IL16", "SPP1", "BTLA", "TNF"),
  "SAunstim_vs_HCunstim" = c("TGFb", "LIGHT"),
  "SAstim_vs_HCstim" = c("TGFb", "IL16", "BTLA", "LIGHT", "BAG")
)

for (comp in comparisons) {
  final_plot_dir <- file.path(base_dir, comp, "Final_Plots")
  setwd(final_plot_dir)
  
  cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = strsplit(comp, "_vs_")[[1]])
  
  pdf(file.path(final_plot_dir, "Signaling_pathway_gene_expression_gdT_cells.pdf"))
  for (pathway in gdT_cell_signaling[[comp]]) {
    plotGeneExpression(cellchat, signaling = pathway, split.by = "datasets", colors.ggplot = TRUE)
  }
  dev.off()
}


