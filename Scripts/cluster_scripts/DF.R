library(data.table)
library(Matrix)
library(Seurat)
library(DoubletFinder)
library(scCustomize)
library(tidyverse)
library(ggplot2)

seurat <- readRDS("/u/scratch/j/jshin/medial_septum_sct/Saves/MS_seurat.rds")
seurat_list <- SplitObject(seurat, split.by = "orig.ident")

seurat_list_df <- lapply(seq_along(seurat_list), function(i) {
  print(i)
  sample <- seurat_list[[i]]
  # (1) Pre-process sample
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample)
  sample <- RunUMAP(sample, dims = 1:30)
  # (2). pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(sample, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  pdf(file = paste0("Plots/QC/doubletfinder/", unique(sample$Sample_GT), ".pdf"),width=3,height=5)
  bcmvn_sample <- find.pK(sweep.stats)
  dev.off()
  mpK <-as.numeric(as.vector(bcmvn_sample$pK[which.max(bcmvn_sample$BCmetric)]))
  # (3). Run DoubletFinder
  DoubletRate = ncol(sample)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
  nExp_poi <- round(DoubletRate*ncol(sample))
  sample <- doubletFinder(sample, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  # (4). Remove doublets
  table(sample@meta.data[,paste0("DF.classifications_0.25_",mpK,"_",nExp_poi)])
  sample$doubletfinder<-sample@meta.data[,paste0("DF.classifications_0.25_",mpK,"_",nExp_poi)]
  
  return(sample)
})

seurat_df <- merge(x = seurat_list_df[[1]], y = seurat_list_df[2:length(seurat_list_df)], merge.data=TRUE)
seurat_df@meta.data[, grepl("pANN", colnames(seurat_df@meta.data))] <- NULL
seurat_df@meta.data[, grepl("DF", colnames(seurat_df@meta.data))] <- NULL
seurat_df@reductions$pca <- seurat@reductions$pca
seurat_df@reductions$umap <- seurat@reductions$umap

save(seurat_df, seurat_list_df, file = "/u/scratch/j/jshin/medial_septum_sct/Saves/MS_DF.rda")





