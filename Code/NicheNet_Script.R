# Loading in the required packages
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
})
setwd('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/NicheNet_Analysis/Plots')

##### Script for NicheNet #####
# Reading in the data
PBMC_filtered_harmony <- readRDS('/Users/tylerjackson/OneDrive - Baylor College of Medicine/Hongjie_Li_Lab_Documents/PBMC_Data_Bin_Su/PBMC_Dataset/Final_Datasets/Final_Annotated_Data/Final_Data_Upload_Analysis/Merged_adata.rds')
table(PBMC_filtered_harmony@meta.data$Final_annotation_broad)
Idents(PBMC_filtered_harmony) <- 'Final_annotation_broad'

# Reading in the NicheNet networks
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

# Define the receiver cell type
receiver = 'gd T-cells'
expressed_genes_receiver <- get_expressed_genes(receiver, PBMC_filtered_harmony, pct = 0.05)

# Define the ligands and receptors for the analysis
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# Define the sender cell types
sender_celltypes <- c("CD4 T-cells", "CD8 T-cells", "Naive B-cells", "NK cells", "Memory/Intermediate B-cells", "Broad cell types", "Treg cells", "Monocytes")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, PBMC_filtered_harmony, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 10233
length(potential_ligands)
## [1] 545
length(potential_ligands_focused)
## [1] 178

# Define the paramters for comparisons between groups, if using pval_adj, no genes are obtained. Using the parameters from the CellChat analysis to be as comparable as possible
condition_oi <-  "SA-stim"
condition_reference <- "HC-stim"

seurat_obj_receiver <- subset(PBMC_filtered_harmony, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "orig.ident",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
head(geneset_oi, 50)

# Create the background gene set from the list of expressed genes in receiver population
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 6501
length(geneset_oi)
## [1] 420 for HC-stim vs HC-unstim
## [1] 333 for SA-stim vs SA-unstim
## [1] 628 for SA-unstim vs HC-unstim
## [1] 805 for SA-stim vs HC-stim

##### Using the sender-focused approach #####
# Calculate ligand activities
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity

# Calculate best upstream ligands
best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

# Start running the sender-focused analysis
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities_all %>% filter(test_ligand %in% potential_ligands_focused)
ligand_activities
best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()
best_upstream_ligands

# Plot the ligand activities on a heatmap <<<< Save this plot
ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% dplyr::select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr

# Target gene plot
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

p_ligand_target

# Receptor plot for the sender-focused approach
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

p_ligand_receptor

# Visualizing the sender cells and their ligands for the sender-focused approach
p_dotplot <- DotPlot(subset(PBMC_filtered_harmony, Final_annotation_broad %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

# Checking the logfoldchange of these ligands in the sender cell types between HC-stim and HC-unstim <<< This also may be a good plot to produce
celltype_order <- levels(Idents(PBMC_filtered_harmony)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
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
                                        legend_title = "LFC")

p_lfc

# Visualizing how the agnostic approach compares to the sender-focused one
#(make_line_plot(ligand_activities = ligand_activities_all,
#                potential_ligands = potential_ligands_focused) +
#    theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))

# Can make an entire summary figure using the code below
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")

# Saving the combined plot below
pdf(paste0(condition_oi, '_vs_', condition_reference, '_Combined_plot_minpct_0.05_logfold_0.25.pdf'), height = 12, width = 20)
combined_plot
dev.off()

# Saving individual figures in the code below
pdf(paste0(condition_oi, '_vs_', condition_reference, '_Separated_plots_minpct_0.05_logfold_0.25.pdf'), height=10, width=12)
p_lfc
p_ligand_aupr
p_ligand_receptor
p_ligand_target
p_dotplot + theme(legend.position = "none",
                  axis.ticks = element_blank(),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(size = 12),
                  axis.text.y = element_text(size = 9),
                  axis.text.x = element_text(size = 9,  angle = 90, hjust = 0))+
  ylab('Expression in Sender')
dev.off()

# Trying a new script out
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))

# Define the ligands and receptors and plot
ligands_oi <- c("TGFB1") # this can be a list of multiple ligands if required
targets_oi <- c("JUN","PMEPA1", "MAP3K5")

active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 4,
                                                      minmax_scaling = TRUE) 


graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")

# To render the graph in RStudio Viewer, uncomment following line of code
DiagrammeR::render_graph(graph_min_max, layout = "tree")

# To export/draw the svg, you need to install DiagrammeRsvg
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))

pdf('Example_interaction_plot_SA-stim_vs_HC-stim.pdf', height = 10, width = 12)
cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))
dev.off()

