#################################################### Loading packages and setting directory for the outputs of NicheNet ##########################################################

### Loading packages and setting directory for the outputs of NicheNet ###
suppressMessages({
  library(Seurat)
  library(SeuratObject)
  library(nichenetr)
  library(SingleCellExperiment)
  library(ggplot2)
  library(scater)
  library(scran)
  library(patchwork)
  library(dplyr)
  library(reticulate)
  library(lattice)
  library(parallel)
  library(cowplot)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(DiagrammeR)
  library(DiagrammeRsvg)
  library(rsvg)
})

# For saving the plots and the results of NicheNet
plot_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/NicheNet_Analysis/Plots/'
RDS_dir <- '/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/NicheNet_Analysis/NicheNet_RDS/'

### Loading in the data and setting the correct ident ###
PBMC_filtered_harmony <- readRDS('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Datasets/Final_Annotated_Data/Final_Data_Upload_Analysis/Merged_adata.rds')
Idents(PBMC_filtered_harmony) <- 'Final_annotation_broad'

#################################################### Read in the NicheNet networks #############################################################

### Read in the NicheNet networks (setting human for this pipeline) ###
organism <- "human"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

#################################################### Configuring sending, receiving cells, and genes used ################################################################

### Configuring the sender and receiver cells and the genes used for the NicheNet pipeline ###
receiver = 'gd T-cells'
expressed_genes_receiver <- get_expressed_genes(receiver, PBMC_filtered_harmony, pct = 0.05)

# Define the ligands and receptors for the analysis
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

# Used in downstream analysis
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# Define the sender cell types
sender_celltypes <- c("CD4 T-cells", "CD8 T-cells", "Naive B-cells", "NK cells", "Memory/Intermediate B-cells", "Broad cell types", "Treg cells", "Monocytes")

# Obtain the expressed genes in each sender cell
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, PBMC_filtered_harmony, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

# Additional checks
length(expressed_genes_sender)
## [1] 10233
length(potential_ligands)
## [1] 545
length(potential_ligands_focused)
## [1] 178

#################################################### DEG analysis function for the NicheNet pipeline ###############################################################

### DEG analysis function for the NicheNet pipeline ###
conditions = list(
  c('HC-stim', 'HC-unstim'),
  c('SA-stim', 'SA-unstim'),
  c('SA-unstim', 'HC-unstim'),
  c('SA-stim', 'HC-stim')
)
print(conditions)

nichenet_pipeline_DEG <- function(data, receiver, conditions = conditions, output_prefix = 'NicheNet') {
  top50_genes <- list() # Storing the results
  geneset_oi_all <- list()
  
  for (i in seq_along(conditions)) {
    condition_oi <- conditions[[i]][1]
    condition_reference <- conditions[[i]][2]
    print(paste0("Reference: ", condition_reference))
    print(paste0("OI: ", condition_oi))
    print(paste0("Processing: ", condition_oi, " vs ", condition_reference))
    
    # Extract receivers
    seurat_obj_receiver <- subset(data, idents = receiver)
    
    # Perform DEG analysis
    DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                      ident.1 = condition_oi, ident.2 = condition_reference,
                                      group.by = "orig.ident",
                                      min.pct = 0.05) %>% rownames_to_column("gene")
    
    # Filter genes
    geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    # Obtain the top 50 genes for the respective comparison
    top50 <- head(geneset_oi, 50)
    
    # Store results in a list
    top50_genes[[paste0(condition_oi, "_vs_", condition_reference)]] <- top50
    
    background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    # Save the results
    output_DEG <- file.path(plot_dir, paste0(output_prefix, "_", condition_oi, "_vs_", condition_reference, ".csv"))
    write.csv(top50, file = output_DEG, row.names = FALSE)
    
    # Store the geneset for downstream analysis
    geneset_oi_all[[paste0(condition_oi, "_vs_", condition_reference)]] <- list(
      "geneset_oi" = geneset_oi,
      "condition_oi" = condition_oi,
      "condition_reference" = condition_reference,
      "background_expressed_genes" = background_expressed_genes
    )
  }
  return(geneset_oi_all) # Return the results
}

geneset_oi_results <- nichenet_pipeline_DEG(data = PBMC_filtered_harmony, receiver = 'gd T-cells', conditions = conditions)

# Save the results as an RDS file
saveRDS(geneset_oi_results, file = file.path(RDS_dir, "geneset_oi.rds"))
geneset_oi_results <- readRDS(file = file.path(RDS_dir, "geneset_oi.rds"))

# Total number of DEGs for each comparison
## [1] 420 for HC-stim vs HC-unstim
## [1] 333 for SA-stim vs SA-unstim
## [1] 628 for SA-unstim vs HC-unstim
## [1] 805 for SA-stim vs HC-stim

#################################################### Function for ligand activities #############################################################

### Function for ligand activities ###
nichenet_ligand_activities <- function(results_list, background_expressed_genes, potential_ligands, ligand_target_matrix, output_dir, output_prefix = "NicheNet_Ligand_Activity") {
  ligand_results <- list()  # Store ligand activities
  
  for (comparison in names(results_list)) {
    # Recall all of the items in the results list
    geneset_oi <- results_list[[comparison]]$geneset_oi
    condition_oi <- results_list[[comparison]]$condition_oi
    condition_reference <- results_list[[comparison]]$condition_reference
    background_expressed_genes <- results_list[[comparison]]$background_expressed_genes
    
    print(paste0("Processing ligand activities for: ", comparison))
    
    # **Calculate Ligand Activities**
    ligand_activities <- predict_ligand_activities(
      geneset = geneset_oi,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands
    ) %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
    
    # **Create Histogram Plot of Ligand Activity**
    p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
      geom_histogram(color="black", fill="darkorange")  + 
      geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
                 color="red", linetype="dashed", size=1) + 
      labs(x="ligand activity (PCC)", y = "# ligands") +
      theme_classic()
    
    # Extract Top Ligands - 100 in this case
    best_upstream_ligands <- ligand_activities %>%
      top_n(100, aupr_corrected) %>%
      arrange(-aupr_corrected) %>%
      pull(test_ligand) %>%
      unique()
    
    # **Save Results**
    output_base <- paste0(output_prefix, "_", condition_oi, "_vs_", condition_reference)
    
    write.csv(ligand_activities, file = file.path(output_dir, paste0(output_base, "_ligand_activities.csv")), row.names = FALSE)
    
    output_pdf <- file.path(output_dir, paste0(output_base, "_ligand_activity_plot",".pdf"))
    ggsave(output_pdf, plot = p_hist_lig_activity, width = 12, height = 8)
    
    # **Store Results**
    ligand_results[[paste0(condition_oi, "_vs_", condition_reference)]] <- list(
      "ligand_activities" = ligand_activities,
      "best_upstream_ligands" = best_upstream_ligands,
      "geneset_oi" = geneset_oi,
      "background_expressed_genes" = background_expressed_genes,
      "condition_oi" = condition_oi,
      "condition_reference" = condition_reference
    )
  }
  return(ligand_results)
}

ligand_activities_results <- nichenet_ligand_activities(results_list = results, background_expressed_genes = background_expressed_genes, 
                                                potential_ligands = potential_ligands, ligand_target_matrix = ligand_target_matrix, output_dir =  plot_dir)

# Save the results as an RDS file
saveRDS(ligand_activities_results, file = file.path(RDS_dir, "ligand_activities_results.rds"))
ligand_activities_results <- readRDS(file = file.path(RDS_dir, "ligand_activities_results.rds"))

#################################################### Ligand-target links function ###########################################################

### Ligand-target links function ###
ligand_target_links <- function(ligand_activities_results, geneset_oi_results, output_dir, output_prefix = "NicheNet_Ligand_Target_Link") {
  ligand_target <- list()
  
  for (comparison in names(geneset_oi_list)) {
    geneset_oi <- geneset_oi_list[[comparison]]$geneset_oi
    condition_oi <- geneset_oi_list[[comparison]]$condition_oi
    condition_reference <- geneset_oi_list[[comparison]]$condition_reference
    background_expressed_genes <- geneset_oi_list[[comparison]]$background_expressed_genes
    
    print(paste0("Processing: ", condition_oi, " vs ", condition_reference))

    for (ligand_comparison in names(ligand_activities_list)) {
      ligand_activities <- ligand_activities_list[[comparison]]$ligand_activities
      best_upstream_ligands <- ligand_activities_list[[comparison]]$best_upstream_ligands
      
      ligand_activities_all <- ligand_activities 
      best_upstream_ligands_all <- best_upstream_ligands
      
      ligand_activities <- ligand_activities_all %>% filter(test_ligand %in% potential_ligands_focused)
      
      best_upstream_ligands <- ligand_activities %>% top_n(100, aupr_corrected) %>% arrange(-aupr_corrected) %>%
        pull(test_ligand) %>% unique()

      ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
        column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
      vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 
      
      p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                           "Prioritized ligands", "Ligand activity", 
                                           legend_title = "AUPR", color = "darkorange") + 
        theme(axis.text.x.top = element_blank(), legend.text = element_text(angle=45, hjust = 1)) + guides(
          fill = guide_colourbar(barwidth =7, barheight = 1)
        )
      active_ligand_target_links_df <- best_upstream_ligands %>%
        lapply(get_weighted_ligand_target_links,
               geneset = geneset_oi,
               ligand_target_matrix = ligand_target_matrix,
               n = 100) %>%
        bind_rows() %>% drop_na()
      
      active_ligand_target_links <- prepare_ligand_target_visualization(
        ligand_target_df = active_ligand_target_links_df,
        ligand_target_matrix = ligand_target_matrix,
        cutoff = 0.33) 
      
      order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
      order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
      
      vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
      
      p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                             color = "purple", legend_title = "Regulatory potential") +
        scale_fill_gradient2(low = "whitesmoke",  high = "purple")
      
      # Saving the plots
      output_base <- paste0(output_prefix, "_", condition_oi, "_vs_", condition_reference)
      
      output_ligand_aupr_pdf <- file.path(output_dir, paste0(output_base, "_ligand_aupr_plot",".pdf"))
      ggsave(output_ligand_aupr_pdf, plot = p_ligand_aupr, width = 12, height = 8)
      
      output_ligand_target_pdf <- file.path(output_dir, paste0(output_base, "_ligand_target_plot",".pdf"))
      ggsave(output_ligand_target_pdf, plot = p_ligand_target, width = 15, height = 10)
      
      # Saving the results
      ligand_target[[paste0(condition_oi, "_vs_", condition_reference)]] <- list(
        "active_ligand_target_links" = active_ligand_target_links,
        "vis_ligand_target" = vis_ligand_target,
        "vis_ligand_aupr" = vis_ligand_aupr,
        "order_ligands" = order_ligands,
        "order_targets" = order_targets,
        "best_upstream_ligands" = best_upstream_ligands,
        "ligand_activities" = ligand_activities,
        "background_expressed_genes" = background_expressed_genes,
        "condition_reference" = condition_reference,
        "condition_oi" = condition_oi
      )
    }
  }
  return(ligand_target)
}

ligand_target_results <- ligand_target_links(ligand_activities_result, geneset_oi_results, plot_dir)

# Save the results as an RDS file
saveRDS(ligand_target_results, file = file.path(RDS_dir, "ligand_targets_results.rds"))
ligand_target_results <- readRDS(file = file.path(RDS_dir, "ligand_targets_results.rds"))

#################################################### Ligand-receptor plots and lfc function ###############################################################
### Ligand-receptor plots and lfc function ###
ligand_receptor_lfc <- function(data, input, output_dir, output_prefix = 'NicheNet_LR_LFC_Plots') {
  Final_NicheNet_Output = list()
  celltype_order <- levels(Idents(data)) 
  
  for (comparison in names(input)) {
    "best_upstream_ligands" = input[[comparison]]$best_upstream_ligands
    "condition_reference" = input[[comparison]]$condition_reference
    "condition_oi" = input[[comparison]]$condition_oi
    "background_expressed_genes" = input[[comparison]]$background_expressed_genes
    "geneset_oi" = input[[comparison]]$geneset_oi
    "active_ligand_target_links" = input[[comparison]]$active_ligand_target_links
    "vis_ligand_target" = input[[comparison]]$vis_ligand_target
    "vis_ligand_aupr" = input[[comparison]]$vis_ligand_aupr
    "order_ligands" = input[[comparison]]$order_ligands
    "order_targets" = input[[comparison]]$order_targets
    "best_upstream_ligands" = input[[comparison]]$best_upstream_ligands
    "ligand_activities" = input[[comparison]]$ligand_activities
    
    
    print(paste0("Processing: ", condition_oi, " vs ", condition_reference))
    
    # Finding ligand-receptor links
    ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
      best_upstream_ligands, expressed_receptors,
      lr_network, weighted_networks$lr_sig) 
    
    # For visualizing network
    vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
      ligand_receptor_links_df,
      best_upstream_ligands,
      order_hclust = "both") 
    
    p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                             y_name = "Ligands", x_name = "Receptors",  
                                             color = "mediumvioletred", legend_title = "Prior interaction potential")
    
    p_dotplot <- DotPlot(subset(data, Final_annotation_broad %in% sender_celltypes),
                         features = rev(best_upstream_ligands), cols = "RdYlBu") + 
      coord_flip() +
      scale_y_discrete(position = "right") + theme(
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 9,  angle = 45, hjust = 0))
    
    # Perform DE of top ligands for visualization
    DE_table_top_ligands <- lapply(
      celltype_order[celltype_order %in% sender_celltypes],
      get_lfc_celltype, 
      seurat_obj = PBMC_filtered_harmony,
      condition_colname = "orig.ident",
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      celltype_col = "Final_annotation_broad",
      min.pct = 0, logfc.threshold = 0,
      features = best_upstream_ligands 
    ) 
    
    DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
      column_to_rownames("gene") 
    
    vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
    
    p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                            "Prioritized ligands", "LFC in Sender",
                                            low_color = "midnightblue", mid_color = "white",
                                            mid = median(vis_ligand_lfc), high_color = "red",
                                            legend_title = "LFC", ) + scale_fill_gradient2(
                                              low = "midnightblue",
                                              mid = "white",
                                              high = "red",
                                              midpoint = median(vis_ligand_lfc),
                                              limits = c(-1,1),
                                              name = "LFC"
                                            )
    p_lfc <- p_lfc + theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),  # Control rotation and alignment
      axis.text.y = element_text(size = 6)  # Adjust y-axis font size
    ) 
    
    # Save these plots as pdf
    output_base <- paste0(output_prefix, "_", condition_oi, "_vs_", condition_reference)
    
    ligand_receptor_pdf <- file.path(output_dir, paste0(output_base, "_ligand_receptor_plot",".pdf"))
    ggsave(ligand_receptor_pdf, plot = p_ligand_receptor, width = 12, height = 8)
    
    ligand_dotplot_pdf <- file.path(output_dir, paste0(output_base, "_ligand_dotplot_plot",".pdf"))
    ggsave(ligand_dotplot_pdf, plot = p_dotplot, width = 15, height = 10)
    
    lfc_pdf <- file.path(output_dir, paste0(output_base, "_ligand_lfc_plot",".pdf"))
    ggsave(lfc_pdf, plot = p_lfc, width = 15, height = 10)
    
    # Save a final object for reference later
    Final_NicheNet_Output[[paste0(condition_oi, "_vs_", condition_reference)]] <- list(
      "vis_ligand_lfc" = vis_ligand_lfc,
      "DE_table_top_ligands" = DE_table_top_ligands,
      "ligand_receptor_links_df" = ligand_receptor_links_df,
      "vis_ligand_receptor_network" = vis_ligand_receptor_network,
      "condition_reference" = condition_reference,
      "condition_oi" = condition_oi,
      "background_expressed_genes" = background_expressed_genes,
      "geneset_oi" = geneset_oi,
      "active_ligand_target_links" = active_ligand_target_links,
      "vis_ligand_target" = vis_ligand_target,
      "vis_ligand_aupr" = vis_ligand_aupr,
      "order_ligands" = order_ligands,
      "order_targets" = order_targets,
      "best_upstream_ligands" = best_upstream_ligands,
      "ligand_activities" = ligand_activities
    )  
  }
  return(Final_NicheNet_Output)
}

Final_NicheNet_Output <- ligand_receptor_lfc(PBMC_filtered_harmony, ligand_target_results, plot_dir)

# Save the Final NicheNet output to RDS file
saveRDS(Final_NicheNet_Output, file = file.path(RDS_dir, 'Final_NicheNet_Output.rds'))
Final_NicheNet_Output <- readRDS(file = file.path(RDS_dir, 'Final_NicheNet_Output.rds'))

#################################################### For customized plots (plotting only specific results) #############################################################
### Plotting the top ligand-target interactions for each comparison ###
HCstim_vs_HCunstim_topLT <- Final_NicheNet_Output[["HC-stim_vs_HC-unstim"]]$vis_ligand_target[c("TNF", "IL1B", "IL10", "TNFSF12"), 1:16, drop = FALSE]
SAstim_vs_SAunstim_topLT <- Final_NicheNet_Output[["SA-stim_vs_SA-unstim"]]$vis_ligand_target[c("IL10", "TGFB1", "TNF"), 1:14, drop = FALSE]
SAunstim_vs_HCunstim_topLT <- Final_NicheNet_Output[["SA-unstim_vs_HC-unstim"]]$vis_ligand_target[c("TGFB1", "TNFSF9", "IL7"), 1:14, drop = FALSE]
SAstim_vs_HCstim_topLT <- Final_NicheNet_Output[["SA-unstim_vs_HC-unstim"]]$vis_ligand_target[c("IL7", "TGFB1", "IL15"), 1:14, drop = FALSE]

# Create a named list of comparisons
comparisons <- list(
  HCstim_vs_HCunstim = HCstim_vs_HCunstim_topLT,
  SAstim_vs_SAunstim = SAstim_vs_SAunstim_topLT,
  SAunstim_vs_HCunstim = SAunstim_vs_HCunstim_topLT,
  SAstim_vs_HCstim = SAstim_vs_HCstim_topLT
)

# Create an empty list to store plots
p_ligand_target <- list()

# Iterate over comparisons and store plots
for (comp_name in names(comparisons)) {
  comp_data <- comparisons[[comp_name]]
  
  p_ligand_target[[comp_name]] <- make_heatmap_ggplot(comp_data, "Prioritized ligands", "Predicted target genes",
                                                      color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke", high = "purple") +  
    theme(
      aspect.ratio = 0.2,                   # Adjust rectangle shape
      axis.title = element_text(size = 16),  # Axis title font size
      axis.text = element_text(size = 14),   # Axis tick labels font size
      legend.title = element_text(size = 14),# Legend title font size
      legend.text = element_text(size = 12)  # Legend tick labels font size
    )
}

# Save the plots to a pdf
for (comp_name in names(comparisons)){
  pdf(file = file.path(plot_dir, paste0(comp_name, '_Separated_plots_minpct_0.05_logfold_0.25_resized_ligandtarget_top.pdf')), height=4, width=10)
  print(p_ligand_target[[comp_name]])
  dev.off()
}

### Plotting the top ligands from the sender cells for each comparison ###
HCstim_vs_HCunstim_topLig <- Final_NicheNet_Output[["HC-stim_vs_HC-unstim"]]$vis_ligand_lfc[c("TNF", "IL1B", "IL10", "TNFSF12"), 1:8, drop = FALSE]
SAstim_vs_SAunstim_topLig <- Final_NicheNet_Output[["SA-stim_vs_SA-unstim"]]$vis_ligand_lfc[c("TNF", "TGFB1", "IL10"), 1:8, drop = FALSE]
SAunstim_vs_HCunstim_topLig <- Final_NicheNet_Output[["SA-unstim_vs_HC-unstim"]]$vis_ligand_lfc[c("IL7", "TNFSF9", "TGFB1"), 1:8, drop = FALSE]
SAstim_vs_HCstim_topLig <- Final_NicheNet_Output[["SA-unstim_vs_HC-unstim"]]$vis_ligand_lfc[c("IL15", "TGFB1", "IL7"), 1:8, drop = FALSE]

# Create a named list of comparisons
comparisons <- list(
  HCstim_vs_HCunstim = HCstim_vs_HCunstim_topLig,
  SAstim_vs_SAunstim = SAstim_vs_SAunstim_topLig,
  SAunstim_vs_HCunstim = SAunstim_vs_HCunstim_topLig,
  SAstim_vs_HCstim = SAstim_vs_HCstim_topLig
)

vis_ligand_lfc <- Final_NicheNet_Output[["HC-stim_vs_HC-unstim"]]$vis_ligand_lfc

# Create an empty list to store plots
p_lfc <- list()

# Iterate over comparisons and store plots
for (comp_name in names(comparisons)) {
  comp_data <- comparisons[[comp_name]]
  
  p_lfc[[comp_name]] <- make_threecolor_heatmap_ggplot(comp_data,
                                                       "Prioritized ligands", "LFC in Sender",
                                                       low_color = "midnightblue", mid_color = "white",
                                                        high_color = "red",
                                                       legend_title = "LFC", ) + scale_fill_gradient2(
                                                         low = "midnightblue",
                                                         mid = "white",
                                                         high = "red",
                                                         limits = c(-0.5,0.5),
                                                         oob = scales::squish, # compress the values to limits of the boundary
                                                         name = "LFC"
                                                       ) + theme(
                                                         aspect.ratio = 0.4,                   # Adjust rectangle shape
                                                         axis.title = element_text(size = 16),  # Axis title font size
                                                         axis.text = element_text(size = 14),   # Axis tick labels font size
                                                         legend.title = element_text(size = 14),# Legend title font size
                                                         legend.text = element_text(size = 12)  # Legend tick labels font size
                                                       )
}

p_lfc[["SAstim_vs_SAunstim"]]

# Save the plots to a pdf
for (comp_name in names(comparisons)){
  pdf(file = file.path(plot_dir, paste0(comp_name, '_Separated_plots_minpct_0.05_logfold_0.25_resized_lfc_top.pdf')), height=7, width=8)
  print(p_lfc[[comp_name]])
  dev.off()
}

### Plotting the interactions between the top ligands and targets for each comparison ###
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))

# Define the ligands and targets for each comparison
ligands_HCstimvsHCunstim <- c('IL10') # Other ligands: TNF, IL1B
targets_HCstimvsHCunstim <- c('STAT1', 'PIM1') 

ligands_SAstimvsSAunstim <- c('TNF') # Other ligands: IL10, TGFB1
targets_SAstimvsSAunstim <- c('PMAIP1', 'TNFAIP3')

ligands_SAunstimvsHCunstim <- c('IL7') # Other ligands: TGFB1, TNFSF9
targets_SAunstimvsHCunstim <- c('IL7R', 'CEBPD') # Other targets: (JUN, ID2), (RORA, ALOX5AP)

ligands_SAstimvsHCstim <- c('IL15') # Other ligands: TGFB1, IL7
targets_SAstimvsHCstim <- c('PIM1', 'TNFAIP3') # Other targets: (JUN, BHLHE40)

ligands <- list(
  HCstim_vs_HCunstim = ligands_HCstimvsHCunstim,
  SAstim_vs_SAunstim = ligands_SAstimvsSAunstim,
  SAunstim_vs_HCunstim = ligands_SAunstimvsHCunstim,
  SAstim_vs_HCstim = ligands_SAstimvsHCstim
)

targets <- list(
  HCstim_vs_HCunstim = targets_HCstimvsHCunstim,
  SAstim_vs_SAunstim = targets_SAstimvsSAunstim,
  SAunstim_vs_HCunstim = targets_SAunstimvsHCunstim,
  SAstim_vs_HCstim = targets_SAstimvsHCstim
)

# Function for plotting ligand-target interaction plots
Map(function(ligands, targets, name) {
  print(paste("Processing", name))
  print(ligands)
  print(targets)
  active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands,
                                                        targets_all = targets,
                                                        weighted_networks = weighted_networks,
                                                        ligand_tf_matrix = ligand_tf_matrix,
                                                        top_n_regulators = 4,
                                                        minmax_scaling = TRUE) 
  
  
  graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                     ligands_all = ligands, targets_all = targets,
                                                     sig_color = "indianred", gr_color = "steelblue") %>% mutate_node_attrs(fontsize = 16, width = 1, fontcolor = '#FFFFFF', fontname = 'Helvetica-Bold', color = 'darkblue', ) 
  
  graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))
  writeLines(graph_svg, "graph.svg")
  rsvg_pdf("graph.svg", file = file.path(plot_dir, paste0(name, ligands, ".pdf")))
}, ligands, targets, names(ligands))

### Creating additional dotplots for genes of interest for each comparison ###
Idents(PBMC_filtered_harmony) <- 'Final_annotation_broad'
gdt_cells <- subset(PBMC_filtered_harmony, subset = Final_annotation_broad == 'gd T-cells')
Idents(gdt_cells) <- 'orig.ident'

# HC-stim vs HC-unstim
gdt_cells_HCstimvsHCunstim <- subset(gdt_cells, subset = orig.ident == 'HC-unstim' | orig.ident == 'HC-stim')
gdt_cells_HCstimvsHCunstim@meta.data$orig.ident <- factor(gdt_cells_HCstimvsHCunstim@meta.data$orig.ident, levels = c("HC-unstim","HC-stim"))
targets_HCstimvsHCunstim <- cbind(genes = colnames(Final_NicheNet_Output[["HC-stim_vs_HC-unstim"]]$vis_ligand_target))
pdf('HCstimvsHCunstim_targets_top10_NicheNet.pdf', height = 4, width = 6)
plot_grid(DotPlot(object = gdt_cells_HCstimvsHCunstim, features = head(targets_HCstimvsHCunstim, 16), group.by = 'orig.ident', scale = FALSE) + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 10), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

# SA-stim vs SA-unstim
gdt_cells_SAstimvsSAunstim <- subset(gdt_cells, subset = orig.ident == 'SA-unstim' | orig.ident == 'SA-stim')
gdt_cells_SAstimvsSAunstim@meta.data$orig.ident <- factor(gdt_cells_SAstimvsSAunstim@meta.data$orig.ident, levels = c("SA-unstim","SA-stim"))
targets_SAstimvsSAunstim <- cbind(genes = colnames(Final_NicheNet_Output[["SA-stim_vs_SA-unstim"]]$vis_ligand_target))
pdf('SAstimvsSAunstim_targets_top10_NicheNet.pdf', height = 4, width = 6)
plot_grid(DotPlot(object = gdt_cells_SAstimvsSAunstim, features = head(targets_SAstimvsSAunstim, 13), group.by = 'orig.ident', scale = FALSE) + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

# SA-unstim vs HC-unstim
gdt_cells_SAunstimvsHCunstim <- subset(gdt_cells, subset = orig.ident == 'HC-unstim' | orig.ident == 'SA-unstim')
gdt_cells_SAunstimvsHCunstim@meta.data$orig.ident <- factor(gdt_cells_SAunstimvsHCunstim@meta.data$orig.ident, levels = c("HC-unstim","SA-unstim"))
targets_SAunstimvsHCunstim <- cbind(genes = colnames(Final_NicheNet_Output[["SA-unstim_vs_HC-unstim"]]$vis_ligand_target))
pdf('SAunstimvsHCunstim_targets_top10_NicheNet.pdf', height = 4, width = 6)
plot_grid(DotPlot(object = gdt_cells_SAunstimvsHCunstim, features = c(head(targets_SAunstimvsHCunstim, 14), 'IL7R'), group.by = 'orig.ident', scale = FALSE) + theme(axis.text.x = element_text(angle = 45, hjust=1, size = 10), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

# SA-stim vs HC-stim
gdt_cells_SAstimvsHCstim <- subset(gdt_cells, subset = orig.ident == 'HC-stim' | orig.ident == 'SA-stim')
gdt_cells_SAstimvsHCstim@meta.data$orig.ident <- factor(gdt_cells_SAstimvsHCstim@meta.data$orig.ident, levels = c("HC-stim","SA-stim"))
targets_SAstimvsHCstim <- cbind(genes = colnames(Final_NicheNet_Output[["SA-stim_vs_HC-stim"]]$vis_ligand_target))
pdf('SAstimvsHCstim_targets_top10_NicheNet.pdf', height = 4, width = 6)
plot_grid(DotPlot(object = gdt_cells_SAstimvsHCstim, features = head(targets_SAstimvsHCstim, 14), group.by = 'orig.ident', scale = FALSE) + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x = element_blank()) + guides(fill=FALSE))
dev.off()

############################################ Session Info ################################################
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
#  [1] rsvg_2.6.1                  DiagrammeRsvg_0.1           DiagrammeR_1.0.11           org.Hs.eg.db_3.20.0         AnnotationDbi_1.68.0        lubridate_1.9.4            
#  [7] forcats_1.0.0               stringr_1.5.1               purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
# [13] tidyverse_2.0.0             cowplot_1.1.3               lattice_0.22-6              reticulate_1.40.0           dplyr_1.1.4                 patchwork_1.3.0            
# [19] scran_1.34.0                scater_1.34.0               scuttle_1.16.0              ggplot2_3.5.1               SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
# [25] Biobase_2.66.0              GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0        
# [31] MatrixGenerics_1.18.1       matrixStats_1.5.0           nichenetr_2.2.0             Seurat_4.3.0                SeuratObject_5.0.2          sp_2.2-0                   

# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0   bitops_1.0-9            httr_1.4.7              RColorBrewer_1.1-3      doParallel_1.0.17       tools_4.4.2             sctransform_0.4.1      
#   [8] backports_1.5.0         R6_2.6.0                lazyeval_0.2.2          uwot_0.2.2              GetoptLong_1.0.5        withr_3.0.2             gridExtra_2.3          
#  [15] fdrtool_1.2.18          progressr_0.15.1        cli_3.6.4               spatstat.explore_3.3-4  fastDummies_1.7.5       spatstat.data_3.1-4     randomForest_4.7-1.2   
#  [22] proxy_0.4-27            ggridges_0.5.6          pbapply_1.7-2           foreign_0.8-88          parallelly_1.42.0       limma_3.62.2            RSQLite_2.3.9          
#  [29] rstudioapi_0.17.1       visNetwork_2.1.2        generics_0.1.3          shape_1.4.6.1           ica_1.0-3               spatstat.random_3.3-2   Matrix_1.7-2           
#  [36] ggbeeswarm_0.7.2        abind_1.4-8             lifecycle_1.0.4         edgeR_4.4.2             recipes_1.1.1           SparseArray_1.6.1       Rtsne_0.17             
#  [43] blob_1.2.4              grid_4.4.2              dqrng_0.4.1             promises_1.3.2          crayon_1.5.3            miniUI_0.1.1.1          beachmat_2.22.0        
#  [50] KEGGREST_1.46.0         metapod_1.14.0          pillar_1.10.1           knitr_1.49              ComplexHeatmap_2.22.0   rjson_0.2.23            future.apply_1.11.3    
#  [57] codetools_0.2-20        glue_1.8.0              V8_6.0.1                spatstat.univar_3.1-1   data.table_1.16.4       vctrs_0.6.5             png_0.1-8              
#  [64] spam_2.11-1             gtable_0.3.6            cachem_1.1.0            gower_1.0.2             xfun_0.50               S4Arrays_1.6.0          mime_0.12              
#  [71] prodlim_2024.06.25      survival_3.8-3          timeDate_4041.110       iterators_1.0.14        hardhat_1.4.1           lava_1.8.1              statmod_1.5.0          
#  [78] bluster_1.16.0          fitdistrplus_1.2-2      ROCR_1.0-11             ipred_0.9-15            nlme_3.1-167            bit64_4.6.0-1           RcppAnnoy_0.0.22       
#  [85] irlba_2.3.5.1           vipor_0.4.7             KernSmooth_2.23-26      rpart_4.1.24            DBI_1.2.3               colorspace_2.1-1        Hmisc_5.2-2            
#  [92] nnet_7.3-20             tidyselect_1.2.1        curl_6.2.0              bit_4.5.0.1             compiler_4.4.2          htmlTable_2.4.3         BiocNeighbors_2.0.1    
#  [99] DelayedArray_0.32.0     plotly_4.10.4           shadowtext_0.1.4        checkmate_2.3.2         scales_1.3.0            caTools_1.18.3          lmtest_0.9-40          
# [106] digest_0.6.37           goftest_1.2-3           spatstat.utils_3.1-2    rmarkdown_2.29          XVector_0.46.0          htmltools_0.5.8.1       pkgconfig_2.0.3        
# [113] base64enc_0.1-3         fastmap_1.2.0           rlang_1.1.5             GlobalOptions_0.1.2     htmlwidgets_1.6.4       UCSC.utils_1.2.0        shiny_1.10.0           
# [120] farver_2.1.2            zoo_1.8-12              jsonlite_1.8.9          BiocParallel_1.40.0     ModelMetrics_1.2.2.2    BiocSingular_1.22.0     magrittr_2.0.3         
# [127] Formula_1.2-5           GenomeInfoDbData_1.2.13 dotCall64_1.2           munsell_0.5.1           Rcpp_1.0.14             viridis_0.6.5           ggnewscale_0.5.0       
# [134] stringi_1.8.4           pROC_1.18.5             zlibbioc_1.52.0         MASS_7.3-64             plyr_1.8.9              listenv_0.9.1           ggrepel_0.9.6          
# [141] deldir_2.0-4            Biostrings_2.74.1       splines_4.4.2           tensor_1.5              hms_1.1.3               circlize_0.4.16         locfit_1.5-9.11        
# [148] igraph_2.1.4            spatstat.geom_3.3-5     RcppHNSW_0.6.0          reshape2_1.4.4          ScaledMatrix_1.14.0     evaluate_1.0.3          renv_1.1.1             
# [155] BiocManager_1.30.25     tzdb_0.4.0              foreach_1.5.2           tweenr_2.0.3            httpuv_1.6.15           RANN_2.6.2              polyclip_1.10-7        
# [162] future_1.34.0           clue_0.3-66             scattermore_1.2         ggforce_0.4.2           rsvd_1.0.5              xtable_1.8-4            e1071_1.7-16           
# [169] RSpectra_0.16-2         later_1.4.1             viridisLite_0.4.2       class_7.3-22            memoise_2.0.1           beeswarm_0.4.0          cluster_2.1.8          
# [176] timechange_0.3.0        globals_0.16.3          caret_7.0-1   







