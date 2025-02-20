##### Analysis of the Gamma-delta T-cell Clusters
###################################################### Importing the necessary libraries and set the working directories ######################################################
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
})

# Set working directories
base_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/'
data_dir <- paste0(base_dir, 'PBMC_Dataset/Final_Datasets/Data_Subsets/')
plot_dir <- paste0(base_dir, 'PBMC_Plots/Correct_PBMC_Analysis/Gdt_Cell_DEG_Analysis/DEG_Analysis_Plots/')
go_plot_dir <- paste0(plot_dir, 'Gdt_Cell_GO_Analysis/Upregulated_Genes_GO/')
go_sheet_dir <- paste0(base_dir, 'GO_Analysis_Sheets/Gdt_Cell_GO/Upregulated_Genes_GO/')

# Load gdt-cell clusters data
Gdt_cell_clusters <- readRDS(file.path(data_dir, '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Datasets/Data_Subsets/Gdt_cells_subset.rds'))

###################################################### Plotting features of the gdT-cells ###################################################### 
# Cell surface markers visualization
pdf(file = file.path(plot_dir, 'TRDC_TRDV_TRGV_expression_gdt_cells.pdf'))
DimPlot(Gdt_cell_clusters, reduction = 'umap', group.by = 'Gdt_cluster_subsets')
FeaturePlot(Gdt_cell_clusters, reduction = 'umap', features = c('TRDC', 'TRDV1', 'TRDV2', 'TRDV3'), order = TRUE)
FeaturePlot(Gdt_cell_clusters, reduction = 'umap', features = c('TRGV2', 'TRGV7', 'TRGV8', 'TRGV10'), order = TRUE)
FeaturePlot(Gdt_cell_clusters, reduction = 'umap', features = c('TRGV3', 'TRGV4', 'TRGV5', 'TRGV9'), order = TRUE)
FeaturePlot(Gdt_cell_clusters, reduction = 'umap', features = c('TRGC1', 'TRGC2', 'TRDC', 'TRGJ2'), order = TRUE)

dotplot_features <- c('TRDC', 'TRDV1', 'TRDV2', 'TRDV3', 'TRGV3', 'TRGV4', 'TRGV5', 'TRGV9')
DotPlot(Gdt_cell_clusters, features = dotplot_features, scale = FALSE, group.by = 'Gdt_cluster_subsets') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  guides(fill = FALSE)
dev.off()

# Analyzing CD27/CD45RA gdT-cell populations in the TRDV2/TRGV9 subcluster
sub_cluster <- 'gdt_cluster_1'
Gdt_cell_clusters_TRDV2_TRGV9 <- subset(Gdt_cell_clusters, subset = Gdt_cluster_subsets == sub_cluster)

# Assign group colors
my_colors <- c('#9FBDD3','#5A6EDC','#F6CC93','#F09924')
Gdt_cell_clusters_TRDV2_TRGV9@meta.data$group.colors <- my_colors[as.integer(Gdt_cell_clusters_TRDV2_TRGV9@meta.data$orig.ident)]

# UMAP plots for TRDV2/TRGV9 cluster
pdf(file = file.path(plot_dir, 'TRGV9TRDV2_gdt_cluster_by_sample_identity.pdf'))
DimPlot(Gdt_cell_clusters_TRDV2_TRGV9, reduction = 'umap', group.by = 'Gdt_cluster_subsets', pt.size = 3)
FeaturePlot(Gdt_cell_clusters_TRDV2_TRGV9, features = 'PTPRC', pt.size = 3) + ggtitle('CD45RA Expression')
FeaturePlot(Gdt_cell_clusters_TRDV2_TRGV9, features = 'CD27', pt.size = 3) + ggtitle('CD27 Expression')
FeaturePlot(Gdt_cell_clusters_TRDV2_TRGV9, reduction = 'umap', features = c('CD27', 'PTPRC'), combine = FALSE, pt.size = 3, blend = TRUE, blend.threshold = 0.5)
DimPlot(Gdt_cell_clusters_TRDV2_TRGV9, group.by = 'orig.ident', cols = my_colors, pt.size = 3)
dev.off()

# Number of gdT-cells per sample group in each subcluster
subclusters <- c('gdt_cluster_1', 'gdt_cluster_2', 'gdt_cluster_3')
gdT_subclusters <- lapply(subclusters, function(cluster) {
  subset(Gdt_cell_clusters, subset = Gdt_cluster_subsets == cluster)
})

names(gdT_subclusters) <- subclusters

lapply(gdT_subclusters, function(cluster) {
  table(cluster@meta.data$orig.ident)
})

pdf(file = file.path(plot_dir, 'IL7R_TGFB1_Expression_Gdt_cells.pdf'))
FeaturePlot(Gdt_cell_clusters, features = c('IL7R'))
FeaturePlot(Gdt_cell_clusters, features = c('TGFB1'))
dev.off()

###################################################### GO and DEG analysis for upregulated genes - comparison of sample groups and gdT-cell subtypes in the gdT-cells  ###################################################### 
# Set cluster identities and ordering
Idents(Gdt_cell_clusters) <- Gdt_cell_clusters@meta.data$orig.ident
Gdt_cell_clusters@meta.data$orig.ident <- factor(
  Gdt_cell_clusters@meta.data$orig.ident, 
  levels = c("HC-unstim", "HC-stim", "SA-unstim", "SA-stim")
)

# Function for DEG and GO analysis between sample groups
perform_deg_go_analysis <- function(obj, ident1, ident2, logfc_threshold = 0.5, ontologies = c("BP", "MF", "CC"), output_prefix) {
  # DEG analysis
  deg_results <- FindMarkers(obj, ident.1 = ident1, ident.2 = ident2, only.pos = TRUE, logfc.threshold = logfc_threshold, test.use = 'wilcox')
  deg_file <- paste0(output_prefix, "_DEG_results.csv")
  write.csv(deg_dir, deg_results, file = deg_file, row.names = TRUE)
  
  # Select top 10 genes based on p-value
  top_genes_10 <- deg_results %>%
    top_n(-10, p_val) %>%
    rownames_to_column("Gene")
  
  
  # Select top 20 genes based on p-value
  top_genes_20 <- deg_results %>%
    top_n(-20, p_val) %>%
    rownames_to_column("Gene")
  
  # Convert gene symbols to ENTREZ IDs
  deg_results <- cbind(Gene = rownames(deg_results), deg_results)
  deg_results_genes <- rownames(deg_results)
  entrez_ids <- bitr(deg_results_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Perform GO analysis for each ontology
  go_results_list <- list()
  for (ont in ontologies) {
      go_results <- enrichGO(
      gene = entrez_ids$ENTREZID,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = ont,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    # Generate GO dotplots
    pdf(paste0(go_plot_dir, output_prefix, "_GO_", ont, ".pdf"))
    print(dotplot(go_results, showCategory = 10, font.size = 12))
    print(dotplot(go_results, showCategory = 20, font.size = 8))
    dev.off()
    
    go_results_list[[ont]] <- go_results
    
    # Save GO results to CSV
    go_file <- paste0(go_sheet_dir, "GO_", ont, "_", output_prefix, ".csv")
    write.csv(data.frame(go_results), file = go_file, row.names = FALSE)
  }
  
  # Plot top 10 genes as a DotPlot
 pdf(paste0(plot_dir, output_prefix, "_top10.pdf"))
 print(DotPlot(
    object = obj, 
    features = top_genes_10$Gene, 
    scale = FALSE, 
    group.by = 'orig.ident'
  ) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
    guides(fill = FALSE))
  dev.off()
  
  # Plot top 20 genes as a DotPlot
  pdf(paste0(plot_dir, output_prefix, "_top20.pdf"))
  print(DotPlot(
    object = obj, 
    features = top_genes_20$Gene, 
    scale = FALSE, 
    group.by = 'orig.ident'
  ) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
    guides(fill = FALSE))
  dev.off()
  
  return(list(deg_results = deg_results, go_results = go_results_list, top_genes_10 = top_genes_10, top_genes_20 = top_genes_20))
}

# Define comparisons for the sample groups in gd T-cells
comparisons_SampleGroups <- list(
  list(ident1 = "HC-stim", ident2 = "HC-unstim", prefix = "HC_stim_vs_HC_unstim_up"),
  list(ident1 = "HC-unstim", ident2 = "HC-stim", prefix = "HC_stim_vs_HC_unstim_down"),
  list(ident1 = "SA-stim", ident2 = "SA-unstim", prefix = "SA_stim_vs_SA_unstim_up"),
  list(ident1 = "SA-unstim", ident2 = "SA-stim", prefix = "SA_stim_vs_SA_unstim_down"),
  list(ident1 = "SA-stim", ident2 = "HC-stim", prefix = "SA_stim_vs_HC_stim_up"),
  list(ident1 = "HC-stim", ident2 = "SA-stim", prefix = "SA_stim_vs_HC_stim_down"),
  list(ident1 = "SA-unstim", ident2 = "HC-unstim", prefix = "SA_unstim_vs_HC_unstim_up"),
  list(ident1 = "HC-unstim", ident2 = "SA-unstim", prefix = "SA_unstim_vs_HC_unstim_down")
)

# Loop through comparisons for the sample groups
Idents(Gdt_cell_clusters) <- 'orig.ident'
results <- lapply(comparisons_SampleGroups, function(comp) {
  perform_deg_go_analysis(
    obj = Gdt_cell_clusters,
    ident1 = comp$ident1,
    ident2 = comp$ident2,
    output_prefix = comp$prefix
  )
})

# Define comparisons for the subtypes of gd T-cells
comparisons_gdtSubsets <- list(
  list(ident1 = "gdt_cluster_1", ident2 = c("gdt_cluster_2", "gdt_cluster_3"), prefix = "gdt_cluster_1_up"),
  list(ident1 = c("gdt_cluster_2", "gdt_cluster_3"), ident2 = "gdt_cluster_1", prefix = "gdt_cluster_1_down"),
  list(ident1 = "gdt_cluster_2", ident2 = c("gdt_cluster_1", "gdt_cluster_3"), prefix = "gdt_cluster_2_up"),
  list(ident1 = c("gdt_cluster_1", "gdt_cluster_3"), ident2 = "gdt_cluster_2", prefix = "gdt_cluster_2_down"),
  list(ident1 = "gdt_cluster_3", ident2 = c("gdt_cluster_1", "gdt_cluster_2"), prefix = "gdt_cluster_3_up"),
  list(ident1 = c("gdt_cluster_1", "gdt_cluster_2"), ident2 = "gdt_cluster_3", prefix = "gdt_cluster_3_down")
)

# Loop through comparisons for the gd T-cell subsets
Idents(Gdt_cell_clusters) <- 'Gdt_cluster_subsets'
results <- lapply(comparisons_gdtSubsets, function(comp) {
  perform_deg_go_analysis(
    obj = Gdt_cell_clusters,
    ident1 = comp$ident1,
    ident2 = comp$ident2,
    output_prefix = comp$prefix
  )
})

# Plotting genes belonging to the top GO terms in each comparison
pdf('Genes_GO_0019221_HC_stim_vs_HC_unstim_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('IRF7', 'NFKBIA', 'CD74', 'SP100', 'PPKACA', 'IFITM1', 'CXCR3', 'ISG15', 'STAT1', 'HSPA1B'), scale = FALSE, group.by = 'orig.ident') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

pdf('Genes_GO_0030030_SA_stim_vs_SA_unstim_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('AKAP9', 'ATP1A1', 'CCDC88C'), scale = FALSE, group.by = 'orig.ident') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

pdf('Genes_GO_0002449_SA_unstim_vs_HC_unstim_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('TGFB1', 'NKG7', 'TRBC2', 'CORO1A', 'GZMB', 'IRF7', 'CD74', 'LYST', 'HLA-DPB1', 'ARID5A', 'HLA-DRB1', 'NCR3'), scale = FALSE, group.by = 'orig.ident') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

pdf('Genes_GO_0002443_SA_stim_vs_HC_stim_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('HLA-A', 'NKG7', 'TGFB1', 'FCGR3A', 'CD2', 'TRBC2', 'LYST', 'TYROBP', 'NCR1', 'ARID5A', 'CD8A', 'GZMB', 'HLA-DRB1', 'HLA-DPB1', 'KLRC4'), scale = FALSE, group.by = 'orig.ident') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

pdf('Genes_GO_0019221_cluster_1_gdt_cluster_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('IL7R','LAPTM5','CXCR5','TNFRSF25','PLCB1','IL18RAP','IFNGR1'), scale = FALSE, group.by = 'Gdt_cluster_subsets') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

pdf('Genes_GO_0002444_cluster_2_gdt_cluster_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('NKG7', 'GZMB', 'FCGR3A', 'LYST', 'TGFB1', 'TYROBP', 'CD8A', 'PRF1', 'HLA-DPB1', 'LAT2', 'NCR1', 'TRBC1', 'KLRC4', 'TBX21', 'HLA-DRB1', 'KLRC2', 'CTSC', 'HLA-DPA1', 'ITGB2', 'ARID5A'), scale = FALSE, group.by = 'Gdt_cluster_subsets') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

pdf('Genes_GO_1903131_cluster_3_gdt_cluster_upregulated.pdf')
plot_grid(DotPlot(object = Gdt_cell_clusters, features =c('LEF1','TCF7','CD27','CAMK4','MYC','RPL22','SOX4','CCR7','TESPA1','CHD7','RHOH','DOCK10','PLCL2','ATM','ST3GAL1','CMTM7','ITM2A'), scale = FALSE, group.by = 'Gdt_cluster_subsets') + 
            theme(axis.text.x = element_text(angle = 45, hjust=1, size = 9), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

# Analysis of other genes in gdT-cells
# Group feature plots
pdf('Gene_expression_analysis.pdf')
features_list <- list(
  c('CD8A', 'CD8B'),
  c('IL7R', 'TGFB1'),
  c('CD3D', 'CD3E')
)

lapply(features_list, function(features) {
  FeaturePlot(Gdt_cell_clusters, features = features, reduction = 'umap', order = TRUE, min.cutoff = 0, max.cutoff = 3)
})

# Make violin plots
lapply(c('CD8A', 'CD8B', 'IL7R', 'TGFB1'), function(feature) {
  VlnPlot(Gdt_cell_clusters, feature, group.by = 'Gdt_cluster_subsets', cols = my_colors)
  VlnPlot(Gdt_cell_clusters, feature, group.by = 'orig.ident', cols = my_colors)
})
dev.off()

#################################### Session Info ####################################
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
#  [1] ggplot2_3.5.1          enrichplot_1.26.6      DOSE_4.0.0             org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0   IRanges_2.40.1         S4Vectors_0.44.0       Biobase_2.66.0        
#  [9] BiocGenerics_0.52.0    clusterProfiler_4.14.4

# loaded via a namespace (and not attached):
#  [1] tidyselect_1.2.1        dplyr_1.1.4             farver_2.1.2            blob_1.2.4              R.utils_2.12.3          Biostrings_2.74.1       lazyeval_0.2.2         
#  [8] fastmap_1.2.0           digest_0.6.37           lifecycle_1.0.4         KEGGREST_1.46.0         tidytree_0.4.6          RSQLite_2.3.9           magrittr_2.0.3         
# [15] compiler_4.4.2          rlang_1.1.5             tools_4.4.2             igraph_2.1.4            data.table_1.16.4       ggtangle_0.0.6          bit_4.5.0.1            
# [22] gson_0.1.0              plyr_1.8.9              RColorBrewer_1.1-3      aplot_0.2.4             BiocParallel_1.40.0     withr_3.0.2             purrr_1.0.4            
# [29] R.oo_1.27.0             grid_4.4.2              GOSemSim_2.32.0         colorspace_2.1-1        GO.db_3.20.0            scales_1.3.0            cli_3.6.4              
# [36] crayon_1.5.3            treeio_1.30.0           generics_0.1.3          rstudioapi_0.17.1       ggtree_3.14.0           httr_1.4.7              reshape2_1.4.4         
# [43] ape_5.8-1               DBI_1.2.3               qvalue_2.38.0           cachem_1.1.0            stringr_1.5.1           zlibbioc_1.52.0         splines_4.4.2          
# [50] parallel_4.4.2          ggplotify_0.1.2         BiocManager_1.30.25     XVector_0.46.0          yulab.utils_0.2.0       vctrs_0.6.5             Matrix_1.7-2           
# [57] jsonlite_1.8.9          gridGraphics_0.5-1      patchwork_1.3.0         bit64_4.6.0-1           ggrepel_0.9.6           tidyr_1.3.1             glue_1.8.0             
# [64] codetools_0.2-20        cowplot_1.1.3           stringi_1.8.4           gtable_0.3.6            GenomeInfoDb_1.42.3     UCSC.utils_1.2.0        munsell_0.5.1          
# [71] tibble_3.2.1            pillar_1.10.1           fgsea_1.32.2            GenomeInfoDbData_1.2.13 R6_2.6.0                lattice_0.22-6          R.methodsS3_1.8.2      
# [78] png_0.1-8               memoise_2.0.1           renv_1.1.1              ggfun_0.1.8             Rcpp_1.0.14             fastmatch_1.1-6         nlme_3.1-167           
# [85] fs_1.6.5                pkgconfig_2.0.3     
