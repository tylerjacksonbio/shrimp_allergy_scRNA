##### Volcano plots for DEGs in the Gamma-delta T-cell clusters
###################################################### Importing the necessary libraries and set the working directories ######################################################
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(ggrepel)
})

# Set the working directories
base_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/'
data_dir <- paste0(base_dir, 'PBMC_Dataset/Final_Datasets/Data_Subsets/')
plot_dir <- paste0(base_dir, 'PBMC_Plots/Correct_PBMC_Analysis/Gdt_Cell_DEG_Analysis/Volcano_Plots/')
deg_dir <- paste0(base_dir, 'PBMC_DEG_Tables/Gdt_Cell_Comparisons/Volcano_Plots_Gene_Lists/')

# Load the Gamma-delta T-cell subclusters RDS
Gdt_cell_clusters <- readRDS(file.path(data_dir, 'Gdt_cells_subset.rds'))
Gdt_cell_clusters

###################################################### Perform the DEG analyses for each comparison and generate volcano plots ###################################################### 
# There will be 8 comparisons total:
# Gdt cluster 2 SA-stim vs SA-unstim
# Compare each of the clusters (3 comparisons total)
# Compare each of the sample groups (4 comparisons total)

# Function for performing the deg analysis
volcano_plot_DEG <- function(
    obj,
    ident1,
    ident2,
    output_prefix = NULL,
    logfc_cutoff = 0.5,
    pval_cutoff = 0.05,
    top_n = 10,
    test_use = "wilcox"
) {
  
  message("Running FindMarkers...")
  
  # Preparing data for plotting on volcano plot
  deg <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = test_use,
    logfc.threshold = -Inf
  )
  
  deg$gene <- rownames(deg)
  
  deg_file <- paste0(output_prefix, "_DEG_results.csv")
  write.csv(deg, paste0(deg_dir,deg_file), row.names = TRUE)
  message("DEG table saved to: ", deg_file)
  
  volcano_df <- deg %>%
    mutate(
      direction = case_when(
        avg_log2FC >= logfc_cutoff & p_val_adj <= pval_cutoff ~ "Upregulated",
        avg_log2FC <= -logfc_cutoff & p_val_adj <= pval_cutoff ~ "Downregulated",
        TRUE ~ "NS"
      ),
      neglog10p = -log10(p_val_adj)
    )
  
  top_up <- volcano_df %>%
    filter(direction == "Upregulated") %>%
    arrange(p_val_adj) %>%
    slice(1:top_n)
  
  top_down <- volcano_df %>%
    filter(direction == "Downregulated") %>%
    arrange(p_val_adj) %>%
    slice(1:top_n)
  
  top_genes <- bind_rows(top_up, top_down)
  
  # Volcano plot code
  p <- ggplot(volcano_df, aes(x = avg_log2FC, y = neglog10p, color = direction)) +
    geom_point(alpha = 0.7, size = 1.5, shape = 16) +
    
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dotted") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dotted") +
    
    scale_color_manual(values = c("Downregulated" = "blue", "NS" = "grey", "Upregulated" = "red")) +
    
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      box.padding = 0.6,
      point.padding = 0.3,
      segment.color = "black",
      max.overlaps = Inf
    ) +
    
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(ident1, " vs ", ident2),
      x = "avg_log2FC",
      y = "-log10(p-value)",
      color = "Legend"
    )
  
  pdf_file <- paste0(output_prefix, "_volcano.pdf")
  ggsave(paste0(plot_dir, pdf_file), p, width = 7, height = 5)
  message("Volcano plot saved to: ", pdf_file)
  
  return(list(
    deg_results = deg,
    volcano_df = volcano_df,
    plot = p
  ))
}

# Perform the comparisons of the sample groups
Idents(Gdt_cell_clusters) <- Gdt_cell_clusters@meta.data$orig.ident

comparisons <- list(
  c("HC-stim", "HC-unstim"),
  c("SA-unstim", "HC-unstim"),
  c("SA-stim", "SA-unstim"),
  c("SA-stim", "HC-stim")
)

for (cmp in comparisons) {
  ident1 <- cmp[1]
  ident2 <- cmp[2]
  
  prefix <- paste0(gsub("-", "", ident1), "_vs_", gsub("-", "", ident2))
  
  volcano_plot_DEG(
    obj = Gdt_cell_clusters,
    ident1 = ident1,
    ident2 = ident2,
    output_prefix = prefix,
    logfc_cutoff = 0.5,
    pval_cutoff = 0.05,
    top_n = 10
  )
}

# Perform the comparisons of each subcluster
Idents(Gdt_cell_clusters) <- Gdt_cell_clusters@meta.data$Gdt_cluster_subsets

comparisons <- list(
  c("gdt_cluster_1", "rest"),
  c("gdt_cluster_2", "rest"),
  c("gdt_cluster_3", "rest"),
)

for (cmp in comparisons) {
  ident1 <- cmp[1]
  ident2 <- cmp[2]
  
  prefix <- paste0(gsub("-", "", ident1), "_vs_", gsub("-", "", ident2))
  
  volcano_plot_DEG(
    obj = Gdt_cell_clusters,
    ident1 = ident1,
    ident2 = NULL,
    output_prefix = prefix,
    logfc_cutoff = 0.5,
    pval_cutoff = 0.05,
    top_n = 10
  )
}

# SA-stim vs SA-unstim comparison specifically for the second gdt subcluster

# Subset the second gdt cell cluster
sub_cluster <- 'gdt_cluster_2'
Gdt_subcluster_2 <- subset(Gdt_cell_clusters, subset = Gdt_cluster_subsets == sub_cluster)

Idents(Gdt_subcluster_2) <- Gdt_subcluster_2@meta.data$orig.ident

# Run the volcano plot code on gdt_subcluster_2
out <- volcano_plot_DEG(
  obj = Gdt_subcluster_2,
  ident1 = "SA-stim",
  ident2 = "SA-unstim",
  output_prefix = "Gdt_subcluster_2_SAstim_vs_SAunstim",
  logfc_cutoff = 0.5,
  pval_cutoff = 0.05,
  top_n = 10
)

#R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.4.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] ggrepel_0.9.6          dplyr_1.1.4            Seurat_5.2.1           SeuratObject_5.0.2     sp_2.2-0               ggplot2_3.5.1          enrichplot_1.26.6      DOSE_4.0.0             org.Hs.eg.db_3.20.0   
[10] AnnotationDbi_1.68.0   IRanges_2.40.1         S4Vectors_0.44.0       Biobase_2.66.0         BiocGenerics_0.52.0    clusterProfiler_4.14.4

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22        splines_4.4.2           later_1.4.1             ggplotify_0.1.2         tibble_3.2.1            R.oo_1.27.0             polyclip_1.10-7         fastDummies_1.7.5       lifecycle_1.0.4        
 [10] globals_0.16.3          lattice_0.22-6          MASS_7.3-64             magrittr_2.0.3          limma_3.62.2            plotly_4.10.4           httpuv_1.6.15           ggtangle_0.0.6          sctransform_0.4.1      
 [19] spam_2.11-1             spatstat.sparse_3.1-0   reticulate_1.40.0       cowplot_1.1.3           pbapply_1.7-2           DBI_1.2.3               RColorBrewer_1.1-3      abind_1.4-8             zlibbioc_1.52.0        
 [28] Rtsne_0.17              purrr_1.0.4             R.utils_2.12.3          yulab.utils_0.2.0       GenomeInfoDbData_1.2.13 irlba_2.3.5.1           listenv_0.9.1           spatstat.utils_3.1-2    tidytree_0.4.6         
 [37] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-2   fitdistrplus_1.2-2      parallelly_1.42.0       codetools_0.2-20        tidyselect_1.2.1        aplot_0.2.4             UCSC.utils_1.2.0       
 [46] farver_2.1.2            matrixStats_1.5.0       spatstat.explore_3.3-4  jsonlite_1.8.9          progressr_0.15.1        ggridges_0.5.6          survival_3.8-3          systemfonts_1.2.1       tools_4.4.2            
 [55] ragg_1.3.3              treeio_1.30.0           ica_1.0-3               Rcpp_1.0.14             glue_1.8.0              gridExtra_2.3           qvalue_2.38.0           GenomeInfoDb_1.42.3     withr_3.0.2            
 [64] BiocManager_1.30.25     fastmap_1.2.0           digest_0.6.37           R6_2.6.0                mime_0.12               gridGraphics_0.5-1      textshaping_1.0.0       colorspace_2.1-1        scattermore_1.2        
 [73] GO.db_3.20.0            tensor_1.5              spatstat.data_3.1-4     RSQLite_2.3.9           R.methodsS3_1.8.2       tidyr_1.3.1             generics_0.1.3          renv_1.1.1              data.table_1.16.4      
 [82] httr_1.4.7              htmlwidgets_1.6.4       uwot_0.2.2              pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4              lmtest_0.9-40           XVector_0.46.0          htmltools_0.5.8.1      
 [91] dotCall64_1.2           fgsea_1.32.2            scales_1.3.0            png_0.1-8               spatstat.univar_3.1-1   ggfun_0.1.8             rstudioapi_0.17.1       reshape2_1.4.4          nlme_3.1-167           
[100] cachem_1.1.0            zoo_1.8-12              stringr_1.5.1           KernSmooth_2.23-26      parallel_4.4.2          miniUI_0.1.1.1          pillar_1.10.1           grid_4.4.2              vctrs_0.6.5            
[109] RANN_2.6.2              promises_1.3.2          xtable_1.8-4            cluster_2.1.8           cli_3.6.4               compiler_4.4.2          rlang_1.1.5             crayon_1.5.3            future.apply_1.11.3    
[118] labeling_0.4.3          plyr_1.8.9              fs_1.6.5                stringi_1.8.4           viridisLite_0.4.2       deldir_2.0-4            BiocParallel_1.40.0     munsell_0.5.1           Biostrings_2.74.1      
[127] lazyeval_0.2.2          spatstat.geom_3.3-5     GOSemSim_2.32.0         Matrix_1.7-2            RcppHNSW_0.6.0          patchwork_1.3.0         bit64_4.6.0-1           future_1.34.0           KEGGREST_1.46.0        
[136] statmod_1.5.0           shiny_1.10.0            ROCR_1.0-11             igraph_2.1.4            memoise_2.0.1           ggtree_3.14.0           fastmatch_1.1-6         bit_4.5.0.1             ape_5.8-1              
[145] gson_0.1.0   

