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
##Typical sex difference gene overlap analysis 
####################################################################

plot_overlap <- function(DEG_df_long, effects = "all", class = "sex") {
  
  # --- 1. Define effects to be used for calculation and plotting ---
  
  # Default list of all possible effects
  all_effects <- c("Typical", "GE", "SCE", "GExSCE", "XCD", "YCD1")
  
  # Determine which effects to use for the rowSums/Shared calculation
  if (length(effects) == 1 && effects == "all") {
    effects_to_use <- all_effects[all_effects != "Typical"]
  } else if (all(effects %in% all_effects)) {
    effects_to_use <- effects
  } else {
    stop("Invalid 'effects' parameter. Must be 'all' or one of: Typical, GE, SCE, GExSCE, XCD, YCD1, YCD2")
  }
  
  # --- 2. Filter initial data based on 'class' parameter ---
  
  if (class == "autosomal") {
    # Filter for non-sex, non-mitochondrial chromosomes (Autosomes)
    DEG_df_filtered <- DEG_df_long %>% 
      dplyr::filter(!CHR %in% c("X", "Y", "MT"))
    
  } else if (class == "sex") {
    # Filter for sex chromosomes (X and Y)
    DEG_df_filtered <- DEG_df_long %>% 
      dplyr::filter(CHR %in% c("X", "Y"))
    
  } else if (class == "all") {
    DEG_df_filtered <- DEG_df_long
  }
  
  # Apply common filters to the already-filtered data
  DEG_df_filtered <- DEG_df_filtered %>%
    dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.1) %>% 
    dplyr::filter(!effect == "YCD2") %>% 
    dplyr::select(gene, logFC, P.Value, adj.P.Val, cell_type, effect) %>% 
    distinct()
  
  # --- 3. Initial Wide Transformation and Filtering ---
  
  sets <- as.data.frame(table(DEG_df_filtered$gene, DEG_df_filtered$effect, DEG_df_filtered$cell_type)) %>% 
    pivot_wider(names_from = Var2, values_from = Freq, values_fill = 0) %>% 
    as.data.frame(check.names = FALSE) %>%
    
    # Ensure all required effect columns are present (fill with 0 if gene/cell_type combo is missing)
    # The columns that exist in 'sets' are filtered against 'all_effects'
    dplyr::select(Var1, Var3, Typical, all_of(effects_to_use)) %>%
    
    # Calculate Typical Overlap filter 
    dplyr::filter(Typical > 0) %>%
    
    # Calculate Shared based on the specified effects_to_use
    dplyr::rowwise() %>% 
    dplyr::mutate(Shared = ifelse(sum(c_across(all_of(effects_to_use))) > 1, 1, 0)) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(Var3) %>%
    dplyr::relocate(Var1, Var3, Typical)
  
  # --- 4. Mutually Exclusive Allocation for Non-Typical Effects ---
  
  # Effects that will be reset to 0 if 'Shared' is 1.
  cols_to_reset <- effects_to_use[effects_to_use %in% colnames(sets)]
  
  for(col in cols_to_reset) {
    sets[[col]] <- ifelse(sets$Shared == 1, 0, sets[[col]])
  }
  
  # --- 5. Summarization and Unique Typical Calculation ---
  
  sets_summary <- sets %>% 
    dplyr::summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(Typical_Unique = Typical - Shared - sum(c_across(all_of(effects_to_use)))) %>%
    dplyr::ungroup()
  
  # --- 6. Proportion Calculation and Reshaping ---
  
  sets_prop <- sets_summary %>%
    # Select the columns used for final plotting/proportion
    # Filter to only include the effects that were actually included in the analysis
    dplyr::select(Var3, Shared, Typical_Unique, all_of(cols_to_reset[cols_to_reset %in% colnames(.)])) %>%
    column_to_rownames(var = "Var3")
  
  # Calculate proportion relative to the Total Typical genes (sum of all displayed categories)
  sets_prop <- t(apply(sets_prop, 1, function(row) row / sum(row))) %>% 
    as.data.frame() %>%
    rownames_to_column(var = "cell_type")
  sets_prop <- left_join(sets_prop, sets_summary[, colnames(sets_summary) %in% c("Var3", "Typical")], 
                         by = c("cell_type" = "Var3")) %>% 
    dplyr::rename(Typical_Total = "Typical") %>% 
    dplyr::mutate(cell_type_label = paste0(cell_type, " (", Typical_Total, ")"))
  
  # Reshape for ggplot
  sets_prop_long <- sets_prop %>% 
    pivot_longer(cols = !c(cell_type, cell_type_label), names_to = "Set", values_to = "Proportion")
  
  # --- 7. Ordering and Final Plotting ---
  
  sets_order <- sets_prop_long %>% 
    dplyr::filter(Set == "Typical_Unique") %>% 
    dplyr::arrange(Proportion) %>% 
    dplyr::pull(cell_type_label)
  
  sets_prop_long$cell_type_label <- factor(sets_prop_long$cell_type_label, levels = sets_order)
  
  # Order the Sets for the stacked bar plot
  set_plot_order <- c("Typical_Unique", "Shared", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2")
  sets_prop_long$Set <- factor(sets_prop_long$Set, 
                               levels = set_plot_order[set_plot_order %in% unique(sets_prop_long$Set)])
  
  # Custom colors (replace with your 'effect_hex_colors' if defined)
  plot_colors <- effect_hex_colors[names(effect_hex_colors) %in% unique(sets_prop_long$Set)]
  plot_colors <- c("Typical_Unique" = "#999999", "Shared" = "#1874CD", plot_colors)
  
  sets_prop_long <- sets_prop_long %>% 
    dplyr::filter(Proportion > 0) %>% 
    dplyr::filter(!is.na(Set))
  
  p1 <- ggplot(sets_prop_long, aes(x = cell_type_label, y = Proportion, fill = Set)) +
    geom_bar(stat = "identity") +
    labs(x = "", y = "Proportion", fill = "Category") +
    ggtitle(paste0(str_to_title(class), " - Effects: ", paste(effects_to_use, collapse = ", "))) + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.95), 
      #plot.title = element_text(size = 35)
      text = element_text(size = 30)
    ) +
    scale_fill_manual(values = plot_colors)
  
  return(p1)
}

p1 <- plot_overlap(DEG_df_long, effect = c("GE", "SCE", "GExSCE"), class = "all")
p2 <- plot_overlap(DEG_df_long, effect = "all", class = "all")
p3 <- plot_overlap(DEG_df_long, effect = c("GE", "SCE", "GExSCE"), class = "sex")
p4 <- plot_overlap(DEG_df_long, effect = "all", class = "sex")
p5 <- plot_overlap(DEG_df_long, effect = c("GE", "SCE", "GExSCE"), class = "autosomal")
p6 <- plot_overlap(DEG_df_long, effect = "all", class = "autosomal")
ggarrange(plotlist = list(p1, p2, p3, p4, p5, p6), nrow = 3, ncol = 2)


sex_effects <- c("GE", "SCE", "GExSCE", "XCD", "YCD1")
overlap_plots <- list()
for(i in 1:length(sex_effects)) {
  
  effect_to_test <- sex_effects[i]
  plots <- plot_overlap(DEG_df_long, effects = effect_to_test, class = "sex") + 
    plot_overlap(DEG_df_long, effects = effect_to_test, class = "autosomal")
  
  overlap_plots[[effect_to_test]] <- plots
}