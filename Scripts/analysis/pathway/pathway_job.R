setwd("/u/project/xyang123/jshin/medial_septum_sct")

.libPaths(c("/u/home/j/jshin/project-xyang123/apps/R/4.3.0", .libPaths))

library(limma)
library(edgeR)
library(data.table)
library(plyr)
library(dplyr)
library(Seurat)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(SingleCellExperiment)
library(magrittr)
library(simplifyEnrichment)
library(stringr)
library(enrichplot)
library(EnhancedVolcano)
library(biomaRt)
library(readxl)
library(UpSetR)
library(tidyr)
library(ggpubr)

data.dir <- "/u/scratch/j/jshin/old-scratch-julianas/medial_septum_sct"
allcells <- c("Astrocyte","Endothelial","Ependymal","Intermediate_Progenitor_Cell","Microglia","Myelinating_Oligodendrocyte","Neuron", "Oligodendrocyte_Progenitor_Cell")

##FUNCTIONS
run_enrich <- function(df, model_to_analyze, sign, currentcell, simplify=TRUE) {
  pathway_results <- list()
  df <- df %>% dplyr::filter(model == model_to_analyze) %>% dplyr::mutate(comparison = factor(comparison))
  
  for(j in 1:nlevels(df$comparison)) {
    currentcomparison <- as.character(levels(df$comparison)[j])

    if(sign == "Up") {
      path_df <- df %>% dplyr::filter(adj.P.Val < 0.05 & logFC > 0 & cell_type == currentcell & comparison == currentcomparison)
    }
    else {
      path_df <- df %>% dplyr::filter(adj.P.Val < 0.05 & logFC < 0 & cell_type == currentcell & comparison == currentcomparison)
    }
    
    if(nrow(path_df) == 0) {next} 
    
    genes <- path_df$gene %>% unlist()
    # entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound = NA)
    # entrezIDs <- as.character(entrezIDs)
    # entrezIDs <- na.omit(entrezIDs)

    go_enrich <- enrichGO(gene = genes,
                      universe = gene_list,
                      OrgDb = org.Mm.eg.db, 
                      keyType = 'SYMBOL',
                      ont = "all",
                      pAdjustMethod = "BH")
    
    #zoe_test <- pathway_enrichment(gene_vector = genes, species = "Mm", GOBP = TRUE, calculate_pathwaywise_medianlog2fc = TRUE, simplify = TRUE, gene_lfc_df = path_df[1:2])
  
    if(nrow(data.frame(go_enrich)) == 0) {next}
    
    go_enrich2 <- clusterProfiler::simplify(go_enrich)
    
    if(nrow(data.frame(go_enrich2)) == 0) {next}
    
    if(simplify == TRUE) {
      results <- data.frame(go_enrich2)
    } else {
      results <- data.frame(go_enrich)
    }
      
    fold <- NULL
    for(r in 1:nrow(results)){
      enriched_genes <- unlist(strsplit(results$geneID[r],"/"))
      enriched_genes <- enriched_genes[enriched_genes %in% path_df$gene]
      fold[r] <- median(path_df[path_df$gene %in% enriched_genes, "logFC"], na.rm = T) }
    
    currentpath <- results
    currentpath$fold <- fold 
    currentpath$sign <- sign
    currentpath$cell_type <- currentcell
    currentpath$comparison <- currentcomparison
    currentpath$model <- model_to_analyze
      
    pathway_results[[model_to_analyze]][[currentcomparison]][[currentcell]] <- currentpath
  }
  return(pathway_results)
}

run_pathway <- function(df, model_to_analyze, cell, simplify=TRUE) {
  if(cell == "All") {
    ct_list <- list()
    for(i in 1:length(allcells)) {
      currentcell <- allcells[i]
      pathway_up_results <- run_enrich(df, model_to_analyze, "Up", currentcell, simplify)
      pathway_down_results <- run_enrich(df, model_to_analyze, "Down", currentcell, simplify)
      ct_list[[currentcell]] <- list(pathway_up_results, pathway_down_results)
    }
    return(ct_list)
    
  } else {
    pathway_up_results <- run_enrich(df, model_to_analyze, "Up", cell)
    pathway_down_results <- run_enrich(df, model_to_analyze, "Down", cell)
    return(list(pathway_up_results, pathway_down_results))
  }
}

bind_cells <- function(x) {
  
  cell_list <- list()
  for(i in 1:length(allcells)) {
    currentcell <- allcells[i]
    df <- do.call("c", unlist(x[[currentcell]], recursive = FALSE, use.names = FALSE)) %>% bind_rows()
    if(nrow(df) == 0) {next}
    cell_list[[currentcell]] <- df[[1]]
    }
  
  pathway_results <- do.call("rbind", cell_list)
  return(pathway_results)
}

load(file = paste0(data.dir, "/Results/Pseudobulk_DEG/DEG_all_models_interaction_results.Rda"))
DEG_df <- do.call(c, unlist(DEG_results, recursive=FALSE)) %>% bind_rows()
DEG_df$model <- factor(DEG_df$model, levels = c("Normative", "FCG", "Dose", "Dose2", "Aneuploidy", "Complete"))
gene_list <- unique(DEG_df$gene)

normative_results <- run_pathway(DEG_df, "Normative", "All", simplify = TRUE)
fcg_results <- run_pathway(DEG_df, "FCG", "All", simplify = TRUE)
dose_results <- run_pathway(DEG_df, "Dose", "All", simplify = TRUE)
dose2_results <- run_pathway(DEG_df, "Dose2", "All", simplify = TRUE)
aneuploidy_results <- run_pathway(DEG_df, "Aneuploidy", "All", simplify = TRUE)
complete_results <- run_pathway(DEG_df, "Complete", "All", simplify = TRUE)

save(normative_results, fcg_results, dose_results, dose2_results, aneuploidy_results, complete_results, file = paste0(data.dir, "/Results/Pseudobulk_DEG/pathway_ALL_interaction2.Rda"))


