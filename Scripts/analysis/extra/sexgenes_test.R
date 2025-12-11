##check sex-specific genes

library(Seurat)
library(dplyr)
library(tidyverse)
library(scCustomize)

seurat <- readRDS("/u/scratch/j/jshin/medial_septum_sct/Saves/For_Figures_new.rds")
seurat_list <- SplitObject(seurat, split.by = "Genotype")

gene_annotation <- read.delim("Gene_annotation/gene_annotation_final.txt")
x_genes <- gene_annotation %>% filter(Chr == "X")
y_genes <- gene_annotation %>% filter(Chr == "Y")

percent_xym <- Percent_Expressing(seurat_object = seurat_list[["XYM"]], features = c("Xist", "Uty"))
percent_xxf <- Percent_Expressing(seurat_object = seurat_list[["XXF"]], features = c("Xist", "Uty"))

percent_x_xym <- Percent_Expressing(seurat_object = seurat_list[["XYM"]], features = x_genes$Symbol)
#percent_x_xym <- percent_x_xym[!rowSums(percent_x_xym) == 0, ]
percent_x_xxf <- Percent_Expressing(seurat_object = seurat_list[["XXF"]], features = x_genes$Symbol)
#percent_x_xxf <- percent_x_xxf[!rowSums(percent_x_xxf) == 0, ]

percent_y_xym <- Percent_Expressing(seurat_object = seurat_list[["XYM"]], features = y_genes$Symbol)
percent_y_xxf <- Percent_Expressing(seurat_object = seurat_list[["XXF"]], features = y_genes$Symbol)
percent_y_xxf %>% summarize(across(everything(), mean))

percents <- data.frame(rbind(percent_x_xym %>% summarize(across(everything(), mean)),
                             percent_x_xxf %>% summarize(across(everything(), mean)),
                             percent_y_xym %>% summarize(across(everything(), mean)),
                             percent_y_xxf %>% summarize(across(everything(), mean))))
rownames(percents) <- c("X_XYM", "X_XXF", "Y_XYM", "Y_XXF")

VlnPlot(seurat_list[["XYM"]], pt.size = 0.1, features = c("Xist", "Uty"), group.by = "cell.type2")
VlnPlot(seurat_list[["XXF"]], pt.size = 0.1, features = c("Xist", "Uty"), group.by = "cell.type2")


#####
seurat_list <- SplitObject(seurat, split.by = "data.sampleName")
samples <- unique(seurat$data.sampleName)

empty_cells <- list()
for(i in 1:length(seurat_list)) {
  
  counts <- seurat_list[[i]]@assays$RNA@counts
  sex_counts <- counts[rownames(counts) %in% c("Xist", "Uty"), ]
  sample <- names(seurat_list[i])
  
  all_columns <- ncol(sex_counts)
  
  # Sum each column and check which columns are empty (i.e., sum equals zero)
  column_sums <- colSums(sex_counts != 0)
  
  # Count the number of columns that are empty (sum equals zero)
  empty_columns <- sum(column_sums == 0)
  
  # Save results 
  results <- data.frame(Sample = sample,
                        Empty_cells = empty_columns,
                        All_cells = all_columns)
  empty_cells[[sample]] <- results
}

empty_cells <- do.call("rbind", empty_cells)

#####

seuratsubset <- subset(seurat, Genotype %in% c("XYM", "XXF"))
metadata <- seuratsubset@meta.data

gene_annotation <- readRDS("/u/project/xyang123/jshin/medial_septum_sct/Data/Derived/gene_annotation.Rds")
gene_annotation2 <- gene_annotation[!duplicated(gene_annotation$mgi_symbol), c(1,3, 7:9)]

sexgenes <- gene_annotation2[gene_annotation2$Class == "Sex-linked", "mgi_symbol"]
counts <- seuratsubset@assays$RNA@counts
sexcounts <- counts[rownames(counts) %in% sexgenes, ]

sex_seurat <- CreateSeuratObject(counts = sexcounts, assay = "RNA", meta.data = metadata)

sex_seurat <- NormalizeData(sex_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sex_seurat <- FindVariableFeatures(sex_seurat, selection.method = "vst")
sex_seurat <- ScaleData(sex_seurat, features = VariableFeatures(object = sex_seurat))
sex_seurat <- RunPCA(sex_seurat, features = VariableFeatures(object = sex_seurat))
sex_seurat <- RunUMAP(sex_seurat, dims = 1:50)
sex_seurat <- FindNeighbors(sex_seurat, dims = 1:50)
sex_seurat <- FindClusters(sex_seurat, resolution = 1)

DimPlot(sex_seurat, reduction = "pca", group.by = "Genotype")
DimPlot(sex_seurat, reduction = "umap", group.by = "Genotype")
DimPlot(sex_seurat, reduction = "umap")




