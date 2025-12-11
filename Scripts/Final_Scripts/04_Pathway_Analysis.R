library(AnnotationHub)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(plyr)
library(dplyr)
library(rrvgo)
library(WebGestaltR)
library(tidyverse)
library(scales)
library(forcats)

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
## Functions to run pathway enrichment on modules
####################################################################
module_pathway_enrichment <- function(DEG_input_df, module_to_test, cutoff) {
  
  module_filter <- modules %>% dplyr::filter(Cluster %in% module_to_test) %>% rownames()
  DEG_df_module <- DEG_input_df %>%
    dplyr::filter(module %in% module_filter) %>% 
    dplyr::mutate(present = 1) %>% 
    dplyr::select(gene, module, present) %>% 
    pivot_wider(names_from = module, values_from = present, values_fill = 0) %>%
    column_to_rownames(var = "gene")
  
  DEG_df_module$sum <- rowSums(DEG_df_module)
  DEG_df_module$keep <- ifelse(DEG_df_module$sum >= (cutoff * (ncol(DEG_df) -1)), "TRUE", "FALSE")
  
  genes <- rownames(DEG_df_module %>% dplyr::filter(keep == TRUE))
  print(length(genes))
  if(length(genes) == 0) {print("Not enough genes - change threshold.")}
  
  pathway_results <- WebGestaltR(enrichMethod = "ORA",
                                 organism = "mmusculus",
                                 enrichDatabase = "geneontology_Biological_Process_noRedundant",
                                 #enrichDatabase = "pathway_REACTOME", 
                                 #interestGeneFolder = "Results/complete_analysis/Pathway/RNA_SD/webgestalt/input/Normative",
                                 interestGene = genes,
                                 interestGeneType = "genesymbol",
                                 referenceSet = "genome",
                                 referenceGeneType = "genesymbol",
                                 #outputDirectory = "Results/complete_analysis/Pathway/RNA_SD/webgestalt/test/",
                                 isOutput = FALSE, 
                                 fdrThr = 0.05)
  
  return(pathway_results)
}

# Function to calculate Jaccard similarity
jaccard_similarity <- function(set1, set2) {
  intersect_len <- length(intersect(set1, set2))
  union_len <- length(union(set1, set2))
  return(intersect_len / union_len)
}

####################################################################
##Run pathway enrichment on modules
####################################################################
modules <- read.delim("Results/complete_analysis/Pathway/metacell_SD/allgene_gene_modules.txt", sep = "\t")

all_pathway_result <- list() 
for(i in 1:length(unique(modules$Cluster))) {
  
  test_module <- i
  
  gene_sets <-  modules %>% dplyr::filter(Cluster == test_module) %>% rownames()
  DEG_df_filter <- DEG_df_long %>% 
    dplyr::filter(adj.P.Val < 0.05 & abs(logFC) >= 0.1) %>%
    dplyr::mutate(module = paste(cell_type, effect, sep = "_")) %>% 
    dplyr::filter(module %in% gene_sets) %>%
    dplyr::select(gene, module) %>% 
    dplyr::mutate(present = 1) %>% 
    pivot_wider(names_from = module, values_from = present, values_fill = 0)
  DEG_df_filter$Module_percent <- rowSums(DEG_df_filter[-1]) / (ncol(DEG_df_filter) - 1)
  
  genes_to_test <- DEG_df_filter %>% dplyr::filter(Module_percent >= 0.20) %>% dplyr::select(gene) %>% unlist(use.names = FALSE) ##originally 0.33 
  
  CP_test <- enrichGO(gene = unique(genes_to_test), OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
  CP_test_simple <- clusterProfiler::simplify(CP_test)
  pathway_result <- CP_test_simple@result
  pathway_result$Cluster <- test_module
  
  all_pathway_result[[i]] <- pathway_result
  
}

filtered_pathway_result <- list()
for(k in 1:length(all_pathway_result)) {
  
  print(k)
  x <- all_pathway_result[[k]]
  pathway_df <- x %>%
    mutate(GeneSet = strsplit(geneID, "/")) 
  if(nrow(x) == 1) {
    filtered_pathway_result[[k]] <- x 
    next
  } 
  
  # Identify and filter similar rows
  to_remove <- c()
  
  for (i in 1:(nrow(pathway_df) - 1)) {
    for (j in (i + 1):nrow(pathway_df)) {
      if (j %in% to_remove) next
      similarity <- jaccard_similarity(pathway_df$GeneSet[[i]], pathway_df$GeneSet[[j]])
      
      if (similarity > 0.5) {  # Adjust threshold as needed
        worse_row <- ifelse(pathway_df$p.adjust[i] < pathway_df$p.adjust[j], j, i) ##originally FDR
        to_remove <- c(to_remove, worse_row)
      }
    }
  }
  
  if(is.null(to_remove)) {
    filtered_pathway_result[[k]] <- x 
    next}
  
  # Keep only unique rows
  filtered_df <- pathway_df[-to_remove, ] %>% 
    dplyr::arrange(p.adjust, desc(Count)) %>% ##originally by FDR and enrichmentRatio %>% 
    dplyr::select(-GeneSet)
  
  filtered_pathway_result[[k]] <- filtered_df
}
filtered_pathway_result <- do.call("rbind", filtered_pathway_result)

filtered_pathway_result_sorted <- filtered_pathway_result %>% 
  dplyr::filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>% 
  slice_min(order_by = p.adjust, n = 50) %>% 
  mutate(present = 1) %>%
  dplyr::select(ID, Description, present, Cluster) %>% 
  pivot_wider(values_from = present, names_from = Cluster, values_fill = 0)
filtered_pathway_result_sorted$Cluster_Count <- rowSums(filtered_pathway_result_sorted[, 3:10]) %>% as.numeric()

# Create a module assignment column to ensure exclusive modules are sorted
filtered_pathway_result_sorted$Cluster_Assignment <- apply(filtered_pathway_result_sorted[ , 3:10], 1, function(x) {
  if(sum(x) == 1) {
    return(which(x == 1))
  } else {
    return(NA)
  }
})

# Sort dataframe: First by number of modules (descending), then by module number (ascending)
filtered_pathway_result_sorted <- filtered_pathway_result_sorted[order(-filtered_pathway_result_sorted$Cluster_Count,
                                                                       filtered_pathway_result_sorted$Cluster_Assignment, na.last = FALSE), ]

write.csv(filtered_pathway_result_sorted, file = "Results/complete_analysis/Pathway/metacell_SD/CP_pathway_results/allgene_filtered_pathway_result_sorted.csv")
save(all_pathway_result, filtered_pathway_result, file = "Results/complete_analysis/Pathway/metacell_SD/CP_pathway_results/allgene_all_pathway_results.rda")
