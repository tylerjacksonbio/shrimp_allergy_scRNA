########################################## This is the Script for Running CellChat ##########################################
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

########################################## Session Info ##########################################
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
# [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     

# other attached packages:
#  [1] reticulate_1.40.0    org.Hs.eg.db_3.20.0  AnnotationDbi_1.68.0 IRanges_2.40.1       S4Vectors_0.44.0     future_1.34.0        lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
# [10] purrr_1.0.4          readr_2.1.5          tidyr_1.3.1          tibble_3.2.1         tidyverse_2.0.0      CellChat_2.1.2       Biobase_2.66.0       BiocGenerics_0.52.0  ggplot2_3.5.1       
# [19] igraph_2.1.4         dplyr_1.1.4          Seurat_5.2.1         SeuratObject_5.0.2   sp_2.2-0            

# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.4.2           later_1.4.1             R.oo_1.27.0             polyclip_1.10-7         ggnetwork_0.5.13        fastDummies_1.7.5      
#   [8] lifecycle_1.0.4         rstatix_0.7.2           doParallel_1.0.17       globals_0.16.3          lattice_0.22-6          MASS_7.3-64             backports_1.5.0        
#  [15] magrittr_2.0.3          plotly_4.10.4           sass_0.4.9              jquerylib_0.1.4         httpuv_1.6.15           NMF_0.28                sctransform_0.4.1      
#  [22] spam_2.11-1             spatstat.sparse_3.1-0   cowplot_1.1.3           pbapply_1.7-2           DBI_1.2.3               RColorBrewer_1.1-3      abind_1.4-8            
#  [29] zlibbioc_1.52.0         Rtsne_0.17              R.utils_2.12.3          yulab.utils_0.2.0       circlize_0.4.16         GenomeInfoDbData_1.2.13 ggrepel_0.9.6          
#  [36] irlba_2.3.5.1           listenv_0.9.1           spatstat.utils_3.1-2    goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-2   fitdistrplus_1.2-2     
#  [43] parallelly_1.42.0       svglite_2.1.3           codetools_0.2-20        DOSE_4.0.0              shape_1.4.6.1           tidyselect_1.2.1        UCSC.utils_1.2.0       
#  [50] farver_2.1.2            matrixStats_1.5.0       spatstat.explore_3.3-4  jsonlite_1.8.9          GetoptLong_1.0.5        BiocNeighbors_2.0.1     Formula_1.2-5          
#  [57] progressr_0.15.1        ggridges_0.5.6          ggalluvial_0.12.5       survival_3.8-3          iterators_1.0.14        systemfonts_1.2.1       foreach_1.5.2          
#  [64] tools_4.4.2             sna_2.8                 ica_1.0-3               Rcpp_1.0.14             glue_1.8.0              gridExtra_2.3           qvalue_2.38.0          
#  [71] GenomeInfoDb_1.42.3     withr_3.0.2             BiocManager_1.30.25     fastmap_1.2.0           digest_0.6.37           timechange_0.3.0        R6_2.6.0               
#  [78] mime_0.12               colorspace_2.1-1        scattermore_1.2         GO.db_3.20.0            tensor_1.5              spatstat.data_3.1-4     RSQLite_2.3.9          
#  [85] R.methodsS3_1.8.2       generics_0.1.3          renv_1.1.1              data.table_1.16.4       FNN_1.1.4.1             httr_1.4.7              htmlwidgets_1.6.4      
#  [92] uwot_0.2.2              pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4              registry_0.5-1          ComplexHeatmap_2.22.0   lmtest_0.9-40          
#  [99] XVector_0.46.0          htmltools_0.5.8.1       carData_3.0-5           dotCall64_1.2           fgsea_1.32.2            clue_0.3-66             scales_1.3.0           
# [106] png_0.1-8               spatstat.univar_3.1-1   rstudioapi_0.17.1       tzdb_0.4.0              rjson_0.2.23            reshape2_1.4.4          coda_0.19-4.1          
# [113] statnet.common_4.11.0   nlme_3.1-167            cachem_1.1.0            zoo_1.8-12              GlobalOptions_0.1.2     KernSmooth_2.23-26      parallel_4.4.2         
# [120] miniUI_0.1.1.1          pillar_1.10.1           grid_4.4.2              vctrs_0.6.5             RANN_2.6.2              ggpubr_0.6.0            promises_1.3.2         
# [127] car_3.1-3               xtable_1.8-4            cluster_2.1.8           cli_3.6.4               compiler_4.4.2          rlang_1.1.5             crayon_1.5.3           
# [134] rngtools_1.5.2          ggsignif_0.6.4          future.apply_1.11.3     plyr_1.8.9              fs_1.6.5                stringi_1.8.4           network_1.19.0         
# [141] viridisLite_0.4.2       deldir_2.0-4            gridBase_0.4-7          BiocParallel_1.40.0     munsell_0.5.1           Biostrings_2.74.1       lazyeval_0.2.2         
# [148] spatstat.geom_3.3-5     GOSemSim_2.32.0         Matrix_1.7-2            RcppHNSW_0.6.0          hms_1.1.3               patchwork_1.3.0         bit64_4.6.0-1          
# [155] KEGGREST_1.46.0         shiny_1.10.0            ROCR_1.0-11             broom_1.0.7             memoise_2.0.1           bslib_0.9.0             fastmatch_1.1-6        
# [162] bit_4.5.0.1       



