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
  
  prefix <- paste0(gsub("-", "_", ident1), "_vs_", gsub("-", "_", ident2))
  
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

