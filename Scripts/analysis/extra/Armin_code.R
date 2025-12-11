library(tidyverse)
library(magrittr)
library(reshape)
library(RColorBrewer)
library(nlme)
library(tibble)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(AnnotationHub)
library(Seurat)
library(org.Mm.eg.db)

########################################
## Prep data 
########################################
seurat <- readRDS(file = "/u/home/j/jshin/scratch/medial_septum_sct/Saves/MS_cellbender_complete_saves/seurat_harmony_sampleid_filtered_annotated.Rds")
sample_metadata <- seurat@meta.data
cell_numbers <- as.data.frame.array(table(seurat$SampleID, seurat$celltype)) %>% rownames_to_column(var = "SampleID")
sample_metadata <- sample_metadata %>% dplyr::select(SampleID, Genotype) %>% 
  left_join(., cell_numbers, by = "SampleID") %>% 
  dplyr::rename(`Inh-IMN` = Inh_IMN,
                `Misc-NT` = Misc_NT)

load(file = "Results/complete_analysis/DEG_processed/DEG_SD_original/DEG_df_RNA_SD.rda")
#add chromosome annotation 
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
ensdb <- ah[["AH73905"]]
columns(ensdb)

gene_info <- select(ensdb,
                    keys = DEG_df$gene,
                    keytype = "SYMBOL",
                    column = c("GENENAME", "ENTREZID", "SEQNAME", "DESCRIPTION"))
gene_info2 <- select(org.Mm.eg.db, 
                     keys = DEG_df$gene,
                     keytype = "SYMBOL", 
                     column = c("ENSEMBL"))
gene_info_all <- left_join(gene_info, gene_info2, by = "SYMBOL") %>% unique()


DEG_df <-DEG_df %>% dplyr::select(gene, logFC, P.Value, adj.P.Val, cell_type, effect)
DEG_df$cell_type <- replace(DEG_df$cell_type, DEG_df$cell_type == "Misc_NT", "Misc-NT")
DEG_df$cell_type <- replace(DEG_df$cell_type, DEG_df$cell_type == "Inh_IMN", "Inh-IMN")
modules <- unique(paste(DEG_df$cell_type, DEG_df$effect, sep = "_"))

ordered_columns <- c("gene", "GENENAME", "ENTREZID", "ENSEMBL", "CHR", "PAR", "DESCRIPTION", 
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
DEG_df_wide <- DEG_df_wide[ordered_columns]

md <- sample_metadata
md$group <- as.factor(md$Genotype) %>%  fct_relevel(., "XXF", "XYF", "XXYF", "XYYF", "XXM", "XYM", "XXYM", "XYYM")
d <- DEG_df_wide

#########################
### set new variables ###
#########################

cell.types <- colnames(md)[c(3:13)]

cell.type.colors <- brewer.pal(n = length(cell.types), name = "Paired")
cell.type.colors <- colorRampPalette(cell.type.colors)(11)
cell.type.colors <- scales::hue_pal()(11)
names(cell.type.colors) <- cell.types

contrasts <- c("Normative", "Gonad", "SexChrom", "Gonad:SexChrom", "addX", "addY1", "addY2")
contrasts.colors <- brewer.pal(n = length(contrasts), name = "Set1")
names(contrasts.colors) <- contrasts

group.colors <- brewer.pal(n = length(levels(md$group)), name = "Set2")
names(group.colors) <- levels(md$group)

##############################################
### DEGs counts across cells and contrasts ###
##############################################

## bar chart fo DEG count per cell_type and contrast combination ##

d %>%
  select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast=1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  arrange(., desc(count)) %>%
  mutate(order=1:nrow(.)) %>%
  mutate(cell_contrast=paste(cell_type, contrast, sep="_")) %>% 
  ggplot(data=., aes(x=reorder(cell_contrast, order), y=count, fill=cell_type)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle=90))

## as a dot plot

d %>%
  select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast = 1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  group_by(cell_type) %>%
  mutate(mean_count = mean(count, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(
    contrast = factor(contrast, levels = c("Normative", "Gonad", "SexChrom", "Xdose", "Ydose")),
    cell_type = reorder(cell_type, mean_count)
  ) %>%
  ggplot(aes(x = contrast, y = cell_type, size = count)) +
  geom_point() +
  scale_size_area(max_size = 10) +
  theme_minimal() +
  labs(x = "Contrast", y = "Cell Type", size = "DEG Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12), axis.text.y = element_text( size=12), axis.title = element_text(size=22), legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16), 
        legend.key.size = unit(1.5, "lines"))


## mixed model to test for DEG differences across cell_type and contrast combinations ##

d %>%
  select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast=1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  lm(data=., count~ contrast + cell_type) %>%
  anova

## is DEG count variation across cell types related to variation in number of cells ? - with and withour granule cells ##

mean_cell_counts <- md %>% 
  select(Astrocyte:Vascular) %>% 
  colMeans %>% 
  as.data.frame %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr:: rename(cell_type=1, mean_count=2) %>%
  arrange(desc(mean_count))

mean_cell_DEG_counts <- d %>%
  select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast = 1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  group_by(cell_type) %>%
  summarise(mean_DEGs = mean(count, na.rm = TRUE)) %>%
  arrange(desc(mean_DEGs))

combined_table <- mean_cell_counts %>%
  inner_join(mean_cell_DEG_counts, by = "cell_type")

combined_table %>% 
  ggplot(data=., aes(x=mean_count, y=mean_DEGs, color=cell_type)) + 
  scale_color_manual(values = cell.type.colors) +
  geom_point(size=7, color="black") +
  geom_point(size=6) +
  xlim(0,100) + # excluding granule cells
  ylim(0, 500) + # excluding granule cells
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text( size=16), axis.title = element_text(size=22), legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16), 
        legend.key.size = unit(1.5, "lines")) +
  labs(x="Mean Cell Count Across Samples", y="Mean DEG Count Across Contrasts")


## plotting proportions of DEGs per chromosome compartment per contrast

d %>%
  filter(CHR!="MT") %>%
  mutate(compartment=fct_other(CHR, keep=c("X", "Y"), other_level="AUT")) %>% 
  select(compartment, matches('adj.P.Val')) %>% 
  mutate(across(where(is.numeric), ~ if_else(. <0.05, 1, 0))) %>% 
  as.data.frame %>% 
  melt(., id.var="compartment") %>%
  mutate(DEG=value) %>%
  mutate(cell.type_contrast=str_remove(variable, "adj.P.Val_")) %>% 
  ggplot(data=., aes(x=cell.type_contrast, y=DEG, fill=compartment)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90))

## testing DEGs for sex chromosome enrichment

data <- d %>%
  dplyr::filter(CHR != "MT") %>%
  dplyr::mutate(compartment = fct_other(CHR, keep = c("X", "Y"), other_level = "AUT")) %>% 
  dplyr::select(compartment, matches('adj.P.Val')) %>%
  mutate(across(where(is.numeric), ~ if_else(. < 0.05, 1, 0))) %>%
  as.data.frame() %>%
  melt(id.var = "compartment") 

# Run chi-squared test per variable
chi_results <- data %>%
  group_by(variable) %>%
  summarise(
    chi_sq_test = list(chisq.test(table(compartment, value))),
    .groups = "drop"
  )

# Extract p-values and test statistics
chi_results <- chi_results %>%
  mutate(
    statistic = sapply(chi_sq_test, function(x) x$statistic),
    p_value = sapply(chi_sq_test, function(x) x$p.value),
    fdr_value = sapply(chi_sq_test, function(x) p.adjust(x$p.value)),
    fdr_value_sig = sapply(chi_sq_test, function(x) (p.adjust(x$p.value)<0.05))
  )

# Display results
chi_results %>%
  dplyr::select(variable, p_value, statistic)

##############################################
### DEGs counts across cells and contrasts ###
##############################################

## bar chart fo DEG count per cell_type and contrast combination ##

d %>%
  dplyr::select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast=1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  arrange(., desc(count)) %>%
  mutate(order=1:nrow(.)) %>%
  mutate(cell_contrast=paste(cell_type, contrast, sep="_")) %>% 
  ggplot(data=., aes(x=reorder(cell_contrast, order), y=count, fill=cell_type)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle=90))

## as a dot plot

d %>%
  dplyr::select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast = 1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  group_by(cell_type) %>%
  mutate(mean_count = mean(count, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(
    contrast = factor(contrast, levels = c("Normative", "Gonad", "SexChrom", "Xdose", "Ydose")),
    cell_type = reorder(cell_type, mean_count)
  ) %>%
  ggplot(aes(x = contrast, y = cell_type, size = count)) +
  geom_point() +
  scale_size_area(max_size = 10) +
  theme_minimal() +
  labs(x = "Contrast", y = "Cell Type", size = "DEG Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12), axis.text.y = element_text( size=12), axis.title = element_text(size=22), legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16), 
        legend.key.size = unit(1.5, "lines"))


## mixed model to test for DEG differences across cell_type and contrast combinations ##

d %>%
  dplyr::select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast=1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  lm(data=., count~ contrast + cell_type) %>%
  anova

## is DEG count variation across cell types related to variation in number of cells ? - with and withour granule cells ##

mean_cell_counts <- md %>% 
  dplyr::select(Astrocyte:Vascular) %>% 
  colMeans %>% 
  as.data.frame %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr:: rename(cell_type=1, mean_count=2) %>%
  arrange(desc(mean_count))

mean_cell_DEG_counts <- d %>%
  dplyr::select(matches('adj.P.Val')) %>%
  apply(MARGIN = 2, function(x) sum(as.numeric(x) < 0.05, na.rm = TRUE)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable") %>%
  dplyr::rename(cell_contrast = 1, count = 2) %>%
  separate(cell_contrast, into = c("Prefix", "cell_type", "contrast"), sep = "_", remove = TRUE) %>%
  group_by(cell_type) %>%
  summarise(mean_DEGs = mean(count, na.rm = TRUE)) %>%
  arrange(desc(mean_DEGs))

combined_table <- mean_cell_counts %>%
  inner_join(mean_cell_DEG_counts, by = "cell_type")

combined_table %>% 
  ggplot(data=., aes(x=mean_count, y=mean_DEGs, color=cell_type)) + 
  scale_color_manual(values = cell.type.colors) +
  geom_point(size=7, color="black") +
  geom_point(size=6) +
  xlim(0,100) + # excluding granule cells
  ylim(0, 500) + # excluding granule cells
  theme(axis.text.x = element_text(size=16), axis.text.y = element_text( size=16), axis.title = element_text(size=22), legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16), 
        legend.key.size = unit(1.5, "lines")) +
  labs(x="Mean Cell Count Across Samples", y="Mean DEG Count Across Contrasts")


## plotting proportions of DEGs per chromosome compartment per contrast

d %>%
  dplyr::filter(CHR!="MT") %>%
  mutate(compartment=fct_other(CHR, keep=c("X", "Y"), other_level="AUT")) %>% 
  dplyr::select(compartment, matches('adj.P.Val')) %>% 
  mutate(across(where(is.numeric), ~ if_else(. <0.05, 1, 0))) %>% 
  as.data.frame %>% 
  melt(., id.var="compartment") %>%
  mutate(DEG=value) %>%
  mutate(cell.type_contrast=str_remove(variable, "adj.P.Val_")) %>% 
  ggplot(data=., aes(x=cell.type_contrast, y=DEG, fill=compartment)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(angle=90))

## testing DEGs for sex chromosome enrichment

data <- d %>%
  dplyr::filter(CHR != "MT") %>%
  mutate(compartment = fct_other(CHR, keep = c("X", "Y"), other_level = "AUT")) %>% 
  dplyr::select(compartment, matches('adj.P.Val')) %>%
  mutate(across(where(is.numeric), ~ if_else(. < 0.05, 1, 0))) %>%
  as.data.frame() %>%
  melt(id.var = "compartment") 

# Run chi-squared test per variable
chi_results <- data %>%
  group_by(variable) %>%
  summarise(
    chi_sq_test = list(chisq.test(table(compartment, value))),
    .groups = "drop"
  )

# Extract p-values and test statistics
chi_results <- chi_results %>%
  mutate(
    statistic = sapply(chi_sq_test, function(x) x$statistic),
    p_value = sapply(chi_sq_test, function(x) x$p.value),
    fdr_value = sapply(chi_sq_test, function(x) p.adjust(x$p.value)),
    fdr_value_sig = sapply(chi_sq_test, function(x) (p.adjust(x$p.value)<0.05))
  )

# Display results
chi_results %>%
  dplyr::select(variable, p_value, statistic)

###############################################################################
### distributions and functional annotations of FC for normative sex effect ###
###############################################################################

# heatmaps for normative sex effects between cell - contrast combinations - stratified for sex chromosomes and autosomes
heatmap.row.column.metadata <- d %>%
  dplyr::select(matches('logFC')) %>% 
  cor(use="pairwise.complete.obs") %>% 
  colnames %>%
    as.data.frame %>% dplyr::rename(row.column.name = 1) %>%
    separate(
      row.column.name, 
      into = c("logFC", "cell.type", "contrast"), 
      sep = "_", 
      extra = "merge" # Keeps "contrast" intact for multi-word contrasts
    ) %>%
  dplyr::filter(contrast=="Normative") %>%
  dplyr::select(-logFC)

# Define color ranges for the heatmap
heatmap_colors <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create annotation for rows and columns
row_annotation <- rowAnnotation(
  CellType = heatmap.row.column.metadata$cell.type,
  Contrast = heatmap.row.column.metadata$contrast,
  col = list(
    CellType = setNames(cell.type.colors, heatmap.row.column.metadata$cell.type),
    Contrast = setNames(rep(contrasts.colors[1], 11), heatmap.row.column.metadata$contrast)
  )
)

col_annotation <- HeatmapAnnotation(
  CellType = heatmap.row.column.metadata$cell.type,
  Contrast = heatmap.row.column.metadata$contrast,
  col = list(
    CellType = setNames(cell.type.colors, heatmap.row.column.metadata$cell.type),
    Contrast = setNames(rep(contrasts.colors[1], 11), heatmap.row.column.metadata$contrast)
  )
)

# Correlation matrix
cor_matrix <- d %>%
  dplyr::filter(CHR %in% c("X", "Y")) %>% # just sex chromosomes
  #dplyr::filter(!(CHR %in% c("X", "Y"))) %>%  # just autosomes
  dplyr::select(matches('logFC') & matches('Normative')) %>%
  cor(use = "pairwise.complete.obs")
#cor_matrix <- cor_matrix[!grepl("Inh_IMN|Misc_NT", rownames(cor_matrix)), ]

# Heatmap
Heatmap(
  cor_matrix,
  name = "Correlation",
  col = heatmap_colors,
  top_annotation = col_annotation,
  left_annotation = row_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)

########################################################################
### exploring patterning of FC across all cell-contrast combinations ###
########################################################################

### quick heatmap of pairwise correlations between all FC vectors for different gene sets

d %>%
  #dplyr::filter(!(CHR %in% c("MT", "X", "Y"))) %>%   # just autosomes
  dplyr::filter(CHR %in% c("X", "Y")) %>% # just sex chromosomes
  #dplyr::filter(CHR %in% c("X", "Y")) %>% # just X chromosomes
  dplyr::select(matches('logFC')) %>% 
  cor(use="pairwise.complete.obs") %>% heatmap


### complex heatmap

# Get row and column metadata 

heatmap.row.column.metadata <- d %>%
  dplyr::select(matches('logFC')) %>% 
  cor(use="pairwise.complete.obs") %>% 
  colnames %>%
  as.data.frame %>% dplyr::rename(row.column.name = 1) %>%
  separate(
    row.column.name, 
    into = c("logFC", "cell.type", "contrast"), 
    sep = "_", 
    extra = "merge" # Keeps "contrast" intact for multi-word contrasts
  ) %>% 
  dplyr::select(-logFC)

# Define color ranges for the heatmap
heatmap_colors <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create annotation for rows and columns
row_annotation <- rowAnnotation(
  CellType = heatmap.row.column.metadata$cell.type,
  Contrast = heatmap.row.column.metadata$contrast,
  col = list(
    CellType = setNames(rep(cell.type.colors, each=7), heatmap.row.column.metadata$cell.type),
    Contrast = setNames(rep(contrasts.colors, 11), heatmap.row.column.metadata$contrast)
  )
)

col_annotation <- HeatmapAnnotation(
  CellType = heatmap.row.column.metadata$cell.type,
  Contrast = heatmap.row.column.metadata$contrast,
  col = list(
    CellType = setNames(rep(cell.type.colors, each=7), heatmap.row.column.metadata$cell.type),
    Contrast = setNames(rep(contrasts.colors, 11), heatmap.row.column.metadata$contrast)
  )
)

# Correlation matrix
cor_matrix <- d %>%
  dplyr::filter(CHR %in% c("X", "Y")) %>%  # just sex chromosomes
  #dplyr::filter(!(CHR %in% c("X", "Y"))) %>%  # just autosomes
  dplyr::select(matches('logFC')) %>%
  cor(use = "pairwise.complete.obs")

# Heatmap
Heatmap(
  cor_matrix,
  name = "Correlation",
  col = heatmap_colors,
  top_annotation = col_annotation,
  left_annotation = row_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)

### for each cell type - which SCT effect best captures the normative FC for autosomal genes ? ###

holder <- matrix(rep(NA, length(cell.types)*(length(contrasts)-1)), nrow=length(cell.types))
rownames(holder) <- cell.types
colnames(holder) <- contrasts[2:7]

for(i in 1:length(cell.types)) {
  
  holder[i,] <- d %>%
    dplyr::filter(CHR %in% c( "X", "Y")) %>%
    #dplyr::filter(!(CHR %in% c("MT", "X", "Y"))) %>%
    dplyr::select(matches('logFC')) %>% 
    dplyr::select(matches(cell.types[i])) %>% 
    cor(use="pairwise.complete.obs") %>%
    extract(1, 2:7)
  
}

# Define color ranges for the heatmap
heatmap_colors <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Create annotation for rows and columns
row_annotation <- rowAnnotation(
  CellType = rownames(holder),
  col = list(
    CellType = setNames(cell.type.colors, rownames(holder))
  )
)

col_annotation <- HeatmapAnnotation(
  Contrast = colnames(holder),
  col = list(
    Contrast = setNames(contrasts.colors[2:7], colnames(holder))
  )
)

# Heatmap
Heatmap(
  holder,
  name = "Correlation",
  col = heatmap_colors,
  top_annotation = col_annotation,
  left_annotation = row_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE
)





