library(data.table)
library(plyr)
library(dplyr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggplot2)
library(ComplexUpset)
library(ggpubr)
library(readxl)
library(tidytext)
library(gtools)
library(GeneOverlap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(corrplot)
library(purrr)
library(EnhancedVolcano)
library(metap)
library(AnnotationHub)
library(magrittr)
library(reshape)
library(nlme)
library(tibble)
library(viridis)
library(emmeans)
library(broom)
library(ggrepel)
library(scales)
library(patchwork)
library(openxlsx)

source("./source/color_palette.R")

####################################################################
##Load DEGs
####################################################################
load(file = "Results/complete_analysis/DEG_processed/DEG_SD_metacell/DEG_df_RNA_SD.rda")

#add chromosome annotation 
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
ensdb <- ah[["AH73905"]]
columns(ensdb)

gene_info <- ensembldb::select(ensdb,
                               keys = DEG_df$gene,
                               keytype = "SYMBOL",
                               column = c("GENENAME", "ENTREZID", "SEQNAME", "DESCRIPTION"))
gene_info2 <- ensembldb::select(org.Mm.eg.db, 
                                keys = DEG_df$gene,
                                keytype = "SYMBOL", 
                                column = c("ENSEMBL"))
gene_info_all <- left_join(gene_info, gene_info2, by = "SYMBOL") %>% unique()


DEG_df <- DEG_df %>% 
  dplyr::select(gene, logFC, P.Value, adj.P.Val, cell_type, effect) %>%
  mutate(cell_type = replace(cell_type, cell_type == "Inh_IMN", "Inh-IMN")) %>%
  mutate(cell_type = replace(cell_type, cell_type == "Misc_NT", "Misc-NT")) %>%
  mutate(cell_type = replace(cell_type, cell_type == "GABA-Chol", "Chol")) %>%
  dplyr::filter(!cell_type == "Misc-NT") %>%
  mutate(module = paste(cell_type, effect, sep = "_"))
modules <- unique(paste(DEG_df$cell_type, DEG_df$effect, sep = "_"))

ordered_columns <- c("gene", "GENENAME", "CHR", 
                     unlist(lapply(modules, function(x) c(paste0("logFC_", x), paste0("P.Value_", x), paste0( "adj.P.Val_", x)))))
par_genes <- c("Mid1", "Erdr1", "Mafl", "Asmtl", "Cd99", "Xg", "Arse", "Sts", "Nlgn4", "Akap17a", "Asmt")

DEG_df_wide <- DEG_df %>%
  pivot_wider(
    names_from = c(cell_type, effect),
    values_from = c(logFC, P.Value, adj.P.Val),
    names_sep = "_"
  )
DEG_df_wide <- DEG_df_wide %>% 
  left_join(gene_info_all, by = c("gene" = "SYMBOL")) %>% 
  dplyr::filter(SEQNAME %in% c(as.character(1:19), "MT", "X", "Y")) %>% 
  dplyr::mutate( PAR = ifelse(gene %in% par_genes, "TRUE", "FALSE")) %>% 
  dplyr::rename("CHR" = "SEQNAME")
DEG_df_wide <- DEG_df_wide[ordered_columns] %>% distinct()

DEG_df_long <- DEG_df %>% 
  left_join(gene_info_all, by = c("gene" = "SYMBOL")) %>% 
  dplyr::filter(SEQNAME %in% c(as.character(1:19), "MT", "X", "Y")) %>% 
  dplyr::mutate( PAR = ifelse(gene %in% par_genes, "TRUE", "FALSE")) %>% 
  dplyr::rename("CHR" = "SEQNAME") %>% 
  dplyr::select(gene, CHR, logFC, P.Value, adj.P.Val, cell_type, effect) %>% distinct()
DEG_df_long$CHR <- factor(DEG_df_long$CHR, levels = c((1:19),"X","Y","MT"))

####################################################################
## Group gene sets based on logFC correlation 
####################################################################
logFC <- DEG_df_long %>% 
  dplyr::mutate(module = paste(cell_type, effect, sep = "_")) %>% 
  dplyr::filter(!is.na(logFC)) %>% 
  dplyr::select(gene, logFC, module, CHR) %>% 
  distinct() %>% 
  pivot_wider(names_from = module, values_from = logFC, values_fill = 0) %>% 
  column_to_rownames(var = "gene")

# Extract annotations from column names
module_names <- colnames(logFC[-1])
annotations <- data.frame(
  celltype = sub("_.*", "", module_names),   # Extract part before the "_"
  effect = sub(".*_", "", module_names)      # Extract part after the "_"
)

# Ensure annotations are factors
annotations$celltype <- factor(annotations$celltype)
annotations$effect <- factor(annotations$effect)

# Create annotations for cell type and effect
names(cell_hex_colors)[names(cell_hex_colors) == "Inh_IMN"] <- "Inh-IMN"
cell_hex_colors <- cell_hex_colors[names(cell_hex_colors) != "Misc_NT"]

rowAnnotation <- rowAnnotation(
  CellType = annotations$celltype,
  Effect = annotations$effect,
  col = list(
    CellType = cell_hex_colors,
    Effect = effect_hex_colors
  )
)

colAnnotation <- HeatmapAnnotation(
  CellType = annotations$celltype,
  Effect = annotations$effect,
  col = list(
    CellType = cell_hex_colors,
    Effect = effect_hex_colors
  )
)

# Create a correlation matrix
cor_matrix <- logFC %>% 
  dplyr::select(-CHR) %>%
  cor(method = "pearson", use = "pairwise.complete.obs") #use removes NA values 

# Plot heatmap
plot <- Heatmap(
  cor_matrix,
  name = "Correlation",
  column_title = "LogFC Correlations (Autosomal)",
  top_annotation = colAnnotation,
  left_annotation = rowAnnotation,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), # Color scale for correlation
  show_row_names = TRUE,
  show_column_names = TRUE
) 

plot <- draw(plot, merge_legends = TRUE)

####################################################################
## Define gene clusters 
####################################################################
hclust_object <- column_dend(plot) %>% as.hclust()
module_groups_h_cut <- cutree(hclust_object, h = 2.6) 
table(module_groups_h_cut)

## 1. Extract and format the Cluster information
# Convert the named vector to a data frame for easier use
cluster_df <- module_groups_h_cut %>%
  as.data.frame() %>%
  dplyr::rename(Cluster = ".") %>%
  dplyr::arrange(Cluster)

# Extract the Cluster vector in the order of the original data columns (cor_matrix)
# Note: module_groups_h_cut is already in the order of the columns in cor_matrix
cluster_annotation_data <- factor(module_groups_h_cut)
names(cluster_annotation_data) <- colnames(cor_matrix) # Ensure names match column names

## 2. Define a color mapping for the Cluster
unique_clusters <- unique(module_groups_h_cut)
cluster_hex_colors <- RColorBrewer::brewer.pal(n = length(unique_clusters), name = "Spectral")
names(cluster_hex_colors) <- unique_clusters

colAnnotation_updated <- HeatmapAnnotation(
  CellType = annotations$celltype,
  Effect = annotations$effect,
  Cluster = cluster_annotation_data,
  col = list(
    CellType = cell_hex_colors,
    Effect = effect_hex_colors,
    Cluster = cluster_hex_colors
  ),
  show_annotation_name = TRUE
)
rowAnnotation_updated <- rowAnnotation(
  CellType = annotations$celltype,
  Effect = annotations$effect,
  Cluster = cluster_annotation_data, # Use the same cluster data
  col = list(
    CellType = cell_hex_colors,
    Effect = effect_hex_colors,
    Cluster = cluster_hex_colors # Use the same color mapping
  )
)

plot_updated <- Heatmap(
  cor_matrix,
  name = "Correlation",
  column_title = "LogFC Correlations (Sex-linked) with Clusters",
  # USE THE UPDATED ANNOTATIONS
  top_annotation = colAnnotation_updated,
  left_annotation = rowAnnotation_updated,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  show_row_names = TRUE,
  show_column_names = TRUE
)

plot_updated <- draw(plot_updated, merge_legends = TRUE)

##save cluster memberships 
write.table(cluster_df, file = "Results/complete_analysis/Pathway/metacell_SD/allgene_gene_modules.txt", quote = FALSE, row.names = TRUE, sep = "\t")