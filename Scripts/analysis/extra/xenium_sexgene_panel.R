##DEG Analysis for Xenium probe

library(Seurat)
library(tidyverse)
library(magrittr)

gene_annotation <- readRDS("/u/project/xyang123/jshin/medial_septum_sct/Data/Derived/gene_annotation.Rds")
gene_annotation2 <- gene_annotation[!duplicated(gene_annotation$mgi_symbol), c(1,3,7:9)]
xenium_panel <- read.csv("/u/project/xyang123/jshin/resources/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv")

load("/u/project/xyang123/jshin/medial_septum_sct/Results/DEG_processed/DEG_v3/DEG_df_celltype2.rda")
MS_DEGs <- DEG_df %>% 
  mutate(module = paste(cell_type, effect, sep = "_"), 
         xenium_panel = ifelse(gene %in% xenium_panel$gene_name, "TRUE", "FALSE")) %>%
  left_join(gene_annotation2, by = c("gene" = "mgi_symbol")) %>%
  filter(Class == "Sex-linked" & adj.P.Val < 0.05 & xenium_panel == "FALSE") %>% 
  select(gene, logFC, module) %>% 
  pivot_wider(names_from = module, values_from = logFC)
MS_DEGs[-1] <- ifelse(is.na(MS_DEGs[-1]), 0, 1)
MS_DEGs$MS_sum <- rowSums(MS_DEGs[ , -1])

load("/u/project/xyang123/jshin/visium_sct/Results/DEG_processed/limma_Spatial_region.Rda")
visium_DEGs <- DEG_df
visium_DEGs$effect <- plyr::mapvalues(visium_DEGs$comparison, 
                                      from = c("Normative", "GonadFemale", "Sex_ChromosomeXX", "Sex_ChromosomeXX:GonadFemale", "additionalX", "additionalY"),
                                      to = c("Normative", "Gonad", "SexChrom", "Gonad:SexChrom", "addX", "addY"))
visium_DEGs <- visium_DEGs %>% 
  mutate(module = paste(region, effect, sep = "_"), 
         xenium_panel = ifelse(gene %in% xenium_panel$gene_name, "TRUE", "FALSE")) %>%
  left_join(gene_annotation2, by = c("gene" = "mgi_symbol")) %>%
  filter(Class == "Sex-linked" & adj.P.Val < 0.05 & xenium_panel == "FALSE") %>% 
  select(gene, logFC, module) %>% 
  pivot_wider(names_from = module, values_from = logFC)
visium_DEGs[-1] <- ifelse(is.na(visium_DEGs[-1]), 0, 1)
visium_DEGs$visium_region_sum <- rowSums(visium_DEGs[ , -1])

load("/u/project/xyang123/jshin/visium_sct/Results/DEG_processed/limma_Spatial_subregion.Rda")
visium_DEGs2 <- DEG_df
visium_DEGs2$effect <- plyr::mapvalues(visium_DEGs2$comparison, 
                                      from = c("Normative", "GonadFemale", "Sex_ChromosomeXX", "Sex_ChromosomeXX:GonadFemale", "additionalX", "additionalY"),
                                      to = c("Normative", "Gonad", "SexChrom", "Gonad:SexChrom", "addX", "addY"))
visium_DEGs2 <- visium_DEGs2 %>% 
  mutate(module = paste(region, effect, sep = "_"), 
         xenium_panel = ifelse(gene %in% xenium_panel$gene_name, "TRUE", "FALSE")) %>%
  left_join(gene_annotation2, by = c("gene" = "mgi_symbol")) %>%
  filter(Class == "Sex-linked" & adj.P.Val < 0.05 & xenium_panel == "FALSE") %>% 
  select(gene, logFC, module) %>% 
  pivot_wider(names_from = module, values_from = logFC)
visium_DEGs2[-1] <- ifelse(is.na(visium_DEGs2[-1]), 0, 1)
visium_DEGs2$visium_subregion_sum <- rowSums(visium_DEGs2[ , -1])

##create full dataframe
all_DEGs <- MS_DEGs %>% 
  full_join(visium_DEGs, by = "gene") %>%
  full_join(visium_DEGs2, by = "gene")
all_DEGs <- all_DEGs %>% 
  select(gene, MS_sum, visium_region_sum, visium_subregion_sum, everything()) %>%            
  arrange(desc(MS_sum), desc(visium_region_sum), desc(visium_subregion_sum)) 
  

