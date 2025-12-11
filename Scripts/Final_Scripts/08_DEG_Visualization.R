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
library(scales)

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
## Control genes dotplots
####################################################################
genesofinterest <- c("Xist", "Tsix", "Eif2s3x", "Ddx3x", "Kdm5c", "Kdm6a", "Kdm5d", "Ddx3y", "Eif2s3y", "Uty")

DEG_df$adj.P.Val <- pmax(DEG_df$adj.P.Val, 1e-200)
DEG_df_filtered <- DEG_df %>% dplyr::filter(gene %in% genesofinterest) %>% mutate(FDR = -log10(adj.P.Val)) %>% 
  dplyr::filter(!is.na(FDR)) %>%
  mutate(effect = factor(effect, levels = c("Normative", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2")),
         cell_type = factor(cell_type, levels = c("GABA", "Glut", "Chol", "Inh-IMN", "Misc-NT", "Astrocyte", "Ependymal", "Oligo", "OPC", "Microglia", "Vascular")),
         gene = factor(gene, levels = genesofinterest))
DEG_df_filtered$logFC2 <- DEG_df_filtered$logFC

DEG_df_filtered %>% 
  ggplot(., aes(y = gene, x = cell_type)) +
  geom_tile(fill = "white") +  xlab("")+ylab("")+ ggtitle("Control genes")+      ## to get the rect filled
  geom_point(aes(colour = logFC2, 
                 size = FDR))  +  ## geom_point for circle illusion
  facet_grid(~effect) + 
  scale_color_gradientn(limits = c(-1.5,1.5), 
                        oob = scales::squish,
                        colours  = colorRampPalette(rev(brewer.pal(11,"RdBu")))(40)) + ## color of the corresponding aes      
  scale_size_continuous(range = c(3, 10)) + ## to tune the size of circles 
  theme_bw()+theme(axis.text.y = element_text(face = "italic", size = 13),
                   axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
                   strip.text.x = element_text(size = 20),
                   plot.margin = unit(c(0.5,1,0.5,1.3), "cm"))+ guides(colour=guide_legend(title="log2FC"),
                                                                       size =guide_legend(title="-log10(FDR)"))


####################################################################
## Count barplots
####################################################################

p1 <- DEG_df_long %>% 
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.1 & effect == "Typical") %>%
  dplyr::mutate(sign = ifelse(logFC > 0, "Up", "Down")) %>%
  {table(.$cell_type, .$effect, .$sign)} %>%
  as.data.frame() %>%
  dplyr::filter(Var2 == "Typical") %>%
  mutate(abs_Freq = Freq, 
         Freq = ifelse(Var3 == "Down", Freq * (-1), Freq)) %>%
  ggplot(., aes(x = reorder_within(Var1, -abs_Freq, Var2), y = Freq, fill = Var3)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~Var2, ncol = 11, scales = "free_x") +
  scale_x_reordered() +
  labs(x = NULL) +
  guides(fill = guide_legend(title = "Direction")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95),
        text = element_text(size = 30))

p2 <- DEG_df_long %>% 
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.1) %>%
  dplyr::mutate(sign = ifelse(logFC > 0, "Up", "Down")) %>%
  {table(.$cell_type, .$effect)} %>%
  as.data.frame() %>%
  dplyr::filter(Var2 == "Typical") %>%
  ggplot(., aes(x = reorder_within(Var1, -Freq, Var2), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~Var2, ncol = 11, scales = "free_x") +
  scale_x_reordered() +
  labs(x = NULL) +
  guides(fill = guide_legend(title = "Cell Type")) +
  scale_fill_manual(values = cell_hex_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95),
        text = element_text(size = 30))

p3 <- DEG_df_long %>% 
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.1) %>%
  dplyr::mutate(sign = ifelse(logFC > 0, "Up", "Down")) %>%
  {table(.$cell_type, .$effect, .$Class)} %>%
  as.data.frame() %>%
  dplyr::filter(Var2 == "Typical") %>%
  ggplot(., aes(x = reorder_within(Var1, -Freq, Var2), y = Freq, fill = Var3)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_wrap(~Var2, ncol = 11, scales = "free_x") +
  scale_x_reordered() +
  labs(x = NULL) +
  guides(fill = guide_legend(title = "Genomic Region")) +
  #scale_fill_manual(values = cell_hex_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.95),
        text = element_text(size = 30))

####################################################################
## Proportion barplots
####################################################################
sc_seurat <- readRDS(file = "/u/scratch/j/jshin/medial_septum_sct/Saves/MS_cellbender_complete_saves/sum_pure_supercell_seurat_harmony_sampleid_filtered_annotated_unfiltered.Rds")

##Proportion counts
cluster.set <- unique(sc_seurat$celltype)
gene.set <- sapply(X = cluster.set, function(c) {
  cells.c <- WhichCells(object = sc_seurat, expression = celltype == c)
  nFeature.c <- sum(rowSums(sc_seurat[['RNA']]@counts[, cells.c ]) != 0) ##number of genes with nonzero counts
  return(nFeature.c)
})

gene.set <- as.data.frame(gene.set)
gene.set$Celltype <- cluster.set
colnames(gene.set)[1] <- "Count"
gene.set$Celltype <- gsub("GABA-Chol", "Chol", gene.set$Celltype)
gene.set$Celltype <- gsub("Inh_IMN", "Inh-IMN", gene.set$Celltype)  

DEG_df %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC) >= 0.1) %>% 
  {table(.$cell_type, .$effect)} %>%
  as.data.frame() %>%
  left_join(gene.set, by = c("Var1" = "Celltype")) %>%
  mutate(Proportion = Freq / Count,
         Var2 = factor(Var2, levels = c("Typical", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2")), 
         Var1 = factor(Var1, levels = c("GABA", "Glut", "Chol", "Inh-IMN", "Astrocyte", "Ependymal", "Oligo", "OPC", "Microglia", "Vascular"))) %>%
  ggplot(., aes(x = reorder_within(Var1, -Proportion, Var2), y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity", position = "identity") +
  facet_grid(~Var2, scales = "free_x") +
  scale_x_reordered() +
  labs(x = NULL) +
  guides(fill = guide_legend(title = "Cell Type")) +
  scale_fill_manual(values = cell_hex_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=.95),
        text = element_text(size = 25))

####################################
## logFC correlation and regression combined
####################################

d <- DEG_df_wide

d$prop.cell.types_normative.DEG <- d %>%
  dplyr::select(contains("Normative") & contains("adj.P.Val")) %>%
  rowwise() %>%
  dplyr::mutate(prop.cell.types_DEG = (sum(c_across(everything()) < 0.05, na.rm = TRUE))/14) %>%
  ungroup() %>%
  pull(prop.cell.types_DEG)

d$mean_normative.log10.adj.P.Val <- d %>%
  dplyr::select(contains("Normative") & contains("adj.P.Val")) %>%
  rowwise() %>%
  dplyr::mutate(mean_normative.adj.P.Val = min(c_across(everything()), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(mean_normative.adj.P.Val = -log10(mean_normative.adj.P.Val)) %>%
  dplyr::mutate(mean_normative.adj.P.Val = replace(mean_normative.adj.P.Val, (mean_normative.adj.P.Val=="Inf" | mean_normative.adj.P.Val>100), 100)) %>%
  pull(mean_normative.adj.P.Val) 

d$mean_normative.logFC <- d %>%
  dplyr::select(contains("Normative") & contains("logFC")) %>%
  rowwise() %>%
  dplyr::mutate(mean_normative.logFC = mean(c_across(everything()), na.rm = TRUE)) %>%
  ungroup() %>%
  pull(mean_normative.logFC)

d$`logFC_Inh-IMN_GExSCE` <- 0
d$`adj.P.Val_Inh-IMN_GExSCE` <- 0

#######################
### AUTOSOMAL GENES ###
#######################
contrasts <- levels(DEG_df$effect)

# CORRELATIONS FOR AUTOSOMAL GENES
sexfc.corrs_aut <- matrix(rep(NA, length(allcells)*(length(contrasts)-1)), nrow=length(allcells))
rownames(sexfc.corrs_aut) <- allcells
colnames(sexfc.corrs_aut) <- contrasts[-1]

for(i in 1:length(allcells)) {
  
  sexfc.corrs_aut[i,] <- d %>%
    #dplyr::filter(CHR %in% c( "X")) %>%  # just X chromosomes
    #dplyr::filter(!gene %in% c("Xist", "Tsix")) %>% # and this removes XIST as well since outlier
    dplyr::filter(!(CHR %in% c("MT", "X", "Y"))) %>% # just autosomes
    dplyr::filter(!prop.cell.types_normative.DEG==0) %>%
    dplyr::select(matches('logFC_')) %>% 
    dplyr::select(matches(allcells[i])) %>% 
    cor(use="pairwise.complete.obs") %>%
    extract(1, -1) 
}


## REGRESSIONS FOR AUTOSOMAL GENES 
sexfc.regr_aut <- matrix(rep(NA, length(allcells)*(length(contrasts)-1)), nrow=length(allcells))
rownames(sexfc.regr_aut) <- allcells
colnames(sexfc.regr_aut) <- contrasts[2:7]

for(i in 1:length(rownames(sexfc.regr_aut))) {
  
  for (j in 1:length(colnames(sexfc.regr_aut))) {
    
    ind.var.to.get <- paste(paste0("logFC_", rownames(sexfc.regr_aut)[i]), colnames(sexfc.regr_aut)[j], sep="_")
    dep.var.to.get <- paste(paste0("logFC_", rownames(sexfc.regr_aut)[i]), "Normative", sep="_")
    
    reg_result <- d %>%
      dplyr::filter(CHR %in% c( "X")) %>%  # just sex chromosomes
      dplyr::filter(!gene %in% c("Xist", "Tsix")) %>% # and this removes XIST as well since outlier
      #dplyr::filter(!(CHR %in% c("MT", "X", "Y"))) %>% # just autosomes
      dplyr:: filter(!prop.cell.types_normative.DEG==0) %>% 
      lm(data=., get(dep.var.to.get) ~ get(ind.var.to.get)) %>% 
      summary %>% 
      use_series("coefficients")
    if(nrow(reg_result) < 2) {
      sexfc.regr_aut[i,j] <- NA
    } else {
      sexfc.regr_aut[i,j] <- reg_result %>% extract(2,1)
    }
  }
}


### COMBINING CORR AND REG
sexfc.corrs_aut_long <- sexfc.corrs_aut  %>%
  as.data.frame %>%
  tibble::rownames_to_column("Cell_Type") %>%
  melt(.,  id.var="Cell_Type") %>% 
  dplyr::rename(Contrast=variable) %>%
  dplyr::rename(Correlation=value)

sexfc.regr_aut_long <- sexfc.regr_aut  %>%
  as.data.frame %>%
  tibble::rownames_to_column("Cell_Type") %>%
  melt(.,  id.var="Cell_Type") %>% 
  dplyr::rename(Contrast=variable) %>%
  dplyr::rename(Regression.slope=value) 

sexfc.corrs.regr_aut_long <- sexfc.corrs_aut_long
sexfc.corrs.regr_aut_long$Regression.slope <- sexfc.regr_aut_long$Regression.slope

### PLOTTING
sexfc.corrs.regr_aut_long_for.plot <- sexfc.corrs.regr_aut_long %>%
  mutate(Cell_Type=fct_relevel(Cell_Type, rev(allcells))) 

ggplot() +
  geom_point(data=sexfc.corrs.regr_aut_long_for.plot, aes(x = Contrast, y = Cell_Type, size = (abs(Correlation)*1.1)), color="black") +
  geom_point(data=sexfc.corrs.regr_aut_long_for.plot, aes(x = Contrast, y = Cell_Type, size = abs(Correlation), color=Regression.slope)) +
  scale_size_area(max_size = 10) +
  labs(x = "Contrast", y = "Cell Type", size = "|correlation|", color="regression slope") +
  scale_color_gradient2(low="blue", mid="white", high="red", limits=c(-1, 1), breaks= c(-1, 0, 1), oob=squish) +
  theme_minimal() +                    
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), 
        axis.text.y = element_text( size=20), 
        axis.title = element_text(size=22), 
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20, margin = margin(b = 5)),
        legend.key.size = unit(1.5, "cm"),
        legend.key.width = unit(1,"cm"))

png("Plots/Final_Figures/Dotplot/metacell_logFC_corr_reg_sexlinked.png", width = 7, height = 8, units = "in", res = 1200)
last_plot()
dev.off()

###############
### X GENES ###
###############
# CORRELATIONS FOR X GENES
sexfc.corrs_X <- matrix(rep(NA, length(allcells)*(length(contrasts)-1)), nrow=length(allcells))
rownames(sexfc.corrs_X ) <- allcells
colnames(sexfc.corrs_X ) <- contrasts[2:7]

for(i in 1:length(allcells)) {
  
  sexfc.corrs_X [i,] <- d %>%
    dplyr::filter(CHR %in% c( "X")) %>%  # just X chromosomes
    dplyr::filter(!gene %in% c("Xist", "Tsix")) %>% # and this removes XIST as well since outlier
    #filter(!(CHR %in% c("MT", "X", "Y"))) %>% # just autosomes
    dplyr::filter(!prop.cell.types_normative.DEG==0) %>%
    dplyr::select(matches('logFC_')) %>% 
    dplyr::select(matches(allcells[i])) %>% 
    cor(use="pairwise.complete.obs") %>%
    extract(1, 2:7)
}

## REGRESSIONS FOR X GENES 
sexfc.regr_X <- matrix(rep(NA, length(allcells)*(length(contrasts)-1)), nrow=length(allcells))
rownames(sexfc.regr_X) <- allcells
colnames(sexfc.regr_X) <- contrasts[2:7]

for(i in 1:length(rownames(sexfc.regr_X))) {
  
  for (j in 1:length(colnames(sexfc.regr_X))) {
    
    ind.var.to.get <- paste(paste0("logFC_", rownames(sexfc.regr_X)[i]), colnames(sexfc.regr_X)[j], sep="_")
    dep.var.to.get <- paste(paste0("logFC_", rownames(sexfc.regr_X)[i]), "Normative", sep="_")
    
    sexfc.regr_X[i,j] <- d %>%
      dplyr::filter(CHR %in% c( "X")) %>%  # just sex chromosomes
      dplyr::filter(!gene %in% c("Xist", "Tisx")) %>% # and this removes XIST as well since outlier
      #filter(!(CHR %in% c("MT", "X", "Y"))) %>% # just autosomes
      dplyr::filter(!prop.cell.types_normative.DEG==0) %>% 
      lm(data=., get(dep.var.to.get) ~ get(ind.var.to.get)) %>% 
      summary %>% 
      use_series("coefficients") %>% 
      extract(2,1)
  }
}


### COMBINING CORR AND REG
sexfc.corrs_X_long <- sexfc.corrs_X  %>%
  as.data.frame %>%
  tibble::rownames_to_column("Cell_Type") %>%
  melt(.,  id.var="Cell_Type") %>% 
  dplyr::rename(Contrast=variable) %>%
  dplyr::rename(Correlation=value)

sexfc.regr_X_long <- sexfc.regr_X  %>%
  as.data.frame %>%
  tibble::rownames_to_column("Cell_Type") %>%
  melt(.,  id.var="Cell_Type") %>% 
  dplyr::rename(Contrast=variable) %>%
  dplyr::rename(Regression.slope=value) 

sexfc.corrs.regr_X_long <- sexfc.corrs_X_long

sexfc.corrs.regr_X_long$Regression.slope <- sexfc.regr_X_long$Regression.slope

sexfc.corrs.regr_X_long_for.plot <- sexfc.corrs.regr_X_long %>%
  mutate(Cell_Type=fct_relevel(Cell_Type, rev(allcells))) 

ggplot() +
  geom_point(data=sexfc.corrs.regr_X_long_for.plot, aes(x = Contrast, y = Cell_Type, size = (abs(Correlation)*1.1)), color="black") +
  geom_point(data=sexfc.corrs.regr_X_long_for.plot, aes(x = Contrast, y = Cell_Type, size = abs(Correlation), color=Regression.slope)) +
  scale_size_area(max_size = 10) +
  labs(x = "Contrast", y = "Cell Type", size = "|correlation|", color="regression slope") +
  scale_color_gradient2(low="blue", mid="white", high="red", limits=c(-1, 1), breaks=c(seq(-1, 0, 1)), oob=squish) +
  theme_minimal() +                    
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20), 
        axis.text.y = element_text( size=20), 
        axis.title = element_text(size=22), 
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20,  margin = margin(b = 5)), 
        legend.key.size = unit(1.5, "lines"))
