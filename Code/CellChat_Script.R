##### Loading in the required packages #####
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

# Define directories and paths
base_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su"
data_dir <- file.path(base_dir, "PBMC_Dataset/Final_Datasets/Final_Annotated_Data")
plot_dir <- file.path(base_dir, "PBMC_Plots/Correct_PBMC_Analysis/Cellchat_Analysis/Comparison_Only_LR_Pairs")
save_dir <- file.path(data_dir, "Cellchat_Analysis_Datasets_Only_LR_Pairs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Load dataset
PBMC_filtered_harmony <- readRDS(file.path(data_dir, "Final_Annotated_Data.rds"))
Idents(PBMC_filtered_harmony) <- "Final_annotation_broad"

# Subset function for CellChat preparation
prepare_cellchat_data <- function(object, condition) {
  subset_obj <- subset(object, subset = orig.ident == condition)
  data_input <- GetAssayData(subset_obj, assay = "RNA", slot = "data")
  meta <- data.frame(group = Idents(subset_obj), row.names = names(Idents(subset_obj)))
  return(list(data = data_input, meta = meta))
}

# Conditions
conditions <- c("HC-unstim", "HC-stim", "SA-unstim", "SA-stim")

# Create CellChat objects
create_and_save_cellchat <- function(condition) {
  prepared_data <- prepare_cellchat_data(PBMC_filtered_harmony, condition)
  cellchat <- createCellChat(object = prepared_data$data, meta = prepared_data$meta, group.by = "group")
  saveRDS(cellchat, file = file.path(save_dir, paste0("cellchat_", condition, "_data.rds")))
  return(cellchat)
}

cellchat_objects <- map(conditions, create_and_save_cellchat)

# Process CellChat objects with a standardized pipeline
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

# Run CellChat processing pipeline
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


# Set global directories
base_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Plots/Correct_PBMC_Analysis/Cellchat_Analysis/Comparison_Only_LR_Pairs"
save_dir <- "/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Datasets/Final_Annotated_Data"

# Set working directory
set_working_dir <- function(sub_dir) {
  dir <- file.path(base_dir, sub_dir)
  setwd(dir)
}

# Save PDF plots
save_pdf <- function(filename, width = 8, height = 10, plot_expr) {
  pdf(filename, width = width, height = height)
  plot_expr()
  dev.off()
}

# Plot gene expressions
plot_gene_expression <- function(cellchat, signaling_list, split_by, output_file) {
  save_pdf(output_file, plot_expr = function() {
    for (signaling in signaling_list) {
      plotGeneExpression(cellchat, signaling = signaling, split.by = split_by, colors.ggplot = TRUE)
    }
  })
}

# Save CellChat objects
save_cellchat_objects <- function(cellchat, comparison_name) {
  saveRDS(cellchat, file.path(save_dir, paste0("cellchat_", comparison_name, ".rds")))
}

# Example: HC-stim vs HC-unstim
comparison_name <- "HCstim_vs_HCunstim"
set_working_dir(comparison_name)
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("HC_unstim", "HC_stim"))

# Plot gene expression
plot_gene_expression(
  cellchat,
  signaling_list = c("MHC-II", "ITGB2", "MHC-I"),
  split_by = "datasets",
  output_file = "Signaling_pathway_gene_expression.pdf"
)

# Save the CellChat object
save_cellchat_objects(cellchat, comparison_name)

# Repeat for other comparisons
comparisons <- list(
  list(name = "SAstim_vs_SAunstim", levels = c("SA_unstim", "SA_stim"), signals = c("MHC-II", "ITGB2", "CLEC", "SELPLG", "SPP1")),
  list(name = "SAunstim_vs_HCunstim", levels = c("HC_unstim", "SA_unstim"), signals = c("PARs", "TGFb", "ITGB2", "MHC-II")),
  list(name = "SAstim_vs_HCstim", levels = c("HC_stim", "SA_stim"), signals = c("TGFb", "BTLA", "TNF", "CLEC"))
)

for (comp in comparisons) {
  set_working_dir(comp$name)
  cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = comp$levels)
  plot_gene_expression(
    cellchat,
    signaling_list = comp$signals,
    split_by = "datasets",
    output_file = "Signaling_pathway_gene_expression.pdf"
  )
  save_cellchat_objects(cellchat, comp$name)
}


