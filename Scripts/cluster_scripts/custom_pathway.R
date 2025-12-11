##Script to run pathway enrichment

organism = "org.Mm.eg.db"
set.seed(99)

library(tidyverse)
library(magrittr)
library(reshape)
library(RColorBrewer)
library(nlme)
library(tibble)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(emmeans)
library(ggplot2)
library(broom)
library(purrr)
library(ggrepel)
library(scales)
library(WebGestaltR)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(organism, character.only = TRUE)
require(DOSE)
library(rrvgo)
library(data.table)
library(GeneSetCluster)
library(forcats)
library(AnnotationHub)
library(org.Mm.eg.db)

source("/u/project/xyang123/jshin/resources/scUtils/PathwayEnrichment/pathway_enrichment_functions_YZ_v2.R")
setwd("/u/project/xyang123/jshin/medial_septum_sct")

################################
################################
################################
#### Read in Data###############
################################
################################
################################
gene_universe <- readRDS("/u/project/xyang123/jshin/cerebellum_sct/Scripts/analysis/Armin/pathway_results/gene_universe.RDS")

modules <- read.delim("Results/complete_analysis/Pathway/RNA_SD/pathway_modules/sexlinked_modules.txt")

load(file = "Results/complete_analysis/DEG_processed/DEG_SD_original/DEG_df_RNA_SD.rda")
DEG_df <- DEG_df %>% 
  mutate(cell_type = replace(cell_type, cell_type == "Inh_IMN", "Inh-IMN")) %>%
  mutate(cell_type = replace(cell_type, cell_type == "Misc_NT", "Misc-NT")) %>%
  mutate(cell_type = replace(cell_type, cell_type == "GABA-Chol", "Chol")) %>%
  dplyr::filter(!cell_type == "Misc-NT") %>% 
  mutate(cell_type = factor(cell_type, levels = c("GABA", "Glut", "Chol", "Inh-IMN", "Astrocyte", "Ependymal", "Oligo", "OPC", "Microglia", "Vascular")), 
         effect = factor(effect, levels = c("Normative", "Gonad", "SexChrom", "Gonad:SexChrom", "addX", "addY1", "addY2")),
         Gene_sets = paste(cell_type, effect, sep = "_")) %>% 
  left_join(modules, by = c("Gene_sets"))

DEG_df <- DEG_df %>% 
  dplyr::filter(!is.na(Module)) %>% 
  #dplyr::filter(effect == "Normative") %>% 
  #dplyr::mutate(Direction = ifelse(logFC > 0, "F", "M"),
  #              Module = paste(cell_type, Direction, sep = "_")) %>% 
  dplyr::filter(adj.P.Val < 0.05) 

DEG_df$avg_logFC <- DEG_df$logFC # Pathway function only takes in avg_logFC but not avg_log2FC.
DEG_df$p_val_adj <- DEG_df$adj.P.Val # Pathway function only takes in avg_logFC but not avg_log2FC.
DEG_df$GENE <- DEG_df$gene
################################
################################
################################
#### PREP DATA FOR OVA##########
################################
################################
################################

# ah <- AnnotationHub()
# query(ah, c("Mus musculus", "Ensembl", "v97"))
# ensdb <- ah[["AH73905"]]
# columns(ensdb)
# 
# gene_info <- select(ensdb,
#                     keys = DEG_df$gene,
#                     keytype = "SYMBOL",
#                     column = c("GENENAME", "ENTREZID", "SEQNAME", "DESCRIPTION"))
# gene_info2 <- select(org.Mm.eg.db,
#                      keys = DEG_df$gene,
#                      keytype = "SYMBOL",
#                      column = c("ENSEMBL"))
# gene_info_all <- left_join(gene_info, gene_info2, by = "SYMBOL") %>% unique()
# 
# 
# DEG_df <-DEG_df %>% dplyr::select(gene, logFC, P.Value, adj.P.Val, cell_type, effect)
# modules <- unique(paste(DEG_df$cell_type, DEG_df$effect, sep = "_"))
# 
# ordered_columns <- c("gene", "GENENAME", "ENTREZID", "ENSEMBL", "CHR", "PAR", "DESCRIPTION",
#                      unlist(lapply(modules, function(x) c(paste0("logFC_", x), paste0("P.Value_", x), paste0( "adj.P.Val_", x)))))
# par_genes <- c("Mid1", "Erdr1", "Mafl", "Asmtl", "Cd99", "Xg", "Arse", "Sts", "Nlgn4", "Akap17a", "Asmt")
# 
# DEG_df_wide <- DEG_df %>%
#   pivot_wider(
#     names_from = c(cell_type, effect),
#     values_from = c(logFC, P.Value, adj.P.Val),
#     names_sep = "_"
#   )
# DEG_df_wide <- DEG_df_wide %>%
#   left_join(gene_info_all, by = c("gene" = "SYMBOL")) %>%
#   dplyr::filter(SEQNAME %in% c(as.character(1:19), "MT", "X", "Y")) %>%
#   dplyr::mutate( PAR = ifelse(gene %in% par_genes, "TRUE", "FALSE")) %>%
#   dplyr::rename("CHR" = "SEQNAME")
# DEG_df_wide <- DEG_df_wide[ordered_columns]
# 
# d <- DEG_df_wide
# 
# cell_DEGs <- d %>%
#   dplyr::select(-c("ENTREZID", "ENSEMBL", "DESCRIPTION")) %>%
#   distinct() %>%
#   column_to_rownames(var = "GENENAME") %>%
#   dplyr::select(matches('logFC|P.Value')) %>%
#   as.data.frame() %>%
#   dplyr::select(matches('Normative')) %>% ##select specific effect here
#   rename_with(
#     ~ ifelse(str_detect(., "logFC_"),
#              paste0("logFC_", str_extract(., "(?<=logFC_)[^_]+")),
#              ifelse(str_detect(., "P.Value_"),
#                     paste0("P.Adjust_", str_extract(., "(?<=P.Value_)[^_]+")),
#                     .))
#   )
# 
# celltype_DEG_sets <- cell_DEGs %>%
#   rownames_to_column(var = "GENE") %>%
#   pivot_longer(
#     cols = -GENE,  # Keep the "Gene" column as is
#     names_to = c("Measure", "Celltype"),  # Separate the column names into "Measure" and "Celltype"
#     names_sep = "_",  # Split the names at the underscore
#     values_to = "Value"  # Put the values (logFC or P.Adjust) into a new column called "Value"
#   ) %>%
#   pivot_wider(
#     names_from = Measure,  # Convert "Measure" into separate columns
#     values_from = Value    # The values in "Value" column should go into the new columns
#   ) %>%
#   dplyr::filter(P.Adjust < 0.05) %>%
#   dplyr::mutate(Direction = ifelse(logFC > 0, "Female", "Male"),  ##add direction
#                 Module = paste(Celltype, Direction, sep = "_"))
# 
# celltype_DEG_sets$avg_logFC <- celltype_DEG_sets$logFC # Pathway function only takes in avg_logFC but not avg_log2FC.
# celltype_DEG_sets$p_val_adj <- celltype_DEG_sets$P.Adjust # Pathway function only takes in avg_logFC but not avg_log2FC.

################################
################################
################################
#### RUN OVA####################
################################
################################
################################

pathway_df <- makePathwayEnrichmentDf(
  DEG_df = DEG_df,
  resources_path = "/u/project/xyang123/jshin/resources/Pathway/", # Change it to your own pathway resource folder. 
  output_Dir = "Results/complete_analysis/Pathway/RNA_SD/custom_script/sexlinked_modules/", # Your choice of output dir
  convertToHuman = FALSE, # FALSE in default. In this case will use mouse-specific pathway DB; if to convert to human gene symbols, then use human-specific pathway DB?
  FDR_threshold=0.05, # Threshold for DEGs to be used for pathway analysis
  logFC_threshold = 0.1, # Threshold for DEGs to be used for pathway analysis
  addlogFC = TRUE,
  MODULE_column="Module", # Can be NULL if want enrichment to be done on Cell_type column.
  total_genes_our_data = gene_universe, # Total genes detected in our dataset. vector of gene names  
  min_max=NULL,
  pval_correction_method = "BH" # Either "Bonferroni" (default) or "BH"
)

saveRDS(pathway_df, file = "Results/complete_analysis/Pathway/RNA_SD/custom_script/sexlinked_modules/pathway_df.RDS")

