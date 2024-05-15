##Limma DEG Script## 

library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(SingleCellExperiment)
library(readxl)
library(ggpubr)
library(ggplotify)
library(scuttle)
library(str2str)
library(dplyr)

############################################################
##Step 1 : Load seurat object and prep metadata 
############################################################

data.dir <- "/u/scratch/j/jshin/old-scratch-julianas/medial_septum_sct"
seurat <- readRDS(file.path(data.dir, "Saves/Pseudobulk_sum_sce_to_seurat.RDS")) ##insert own file path here

##Edit metadata 
allcells <- unique(seurat$cell.type)

seurat$Xdose <- plyr::mapvalues(seurat$SexChrom,
                                       from = c("XX", "XXY", "XY", "XYY"),
                                       to = c(2, 2, 1 ,1)) %>% factor(levels = c(1, 2))
seurat$Xdose_cont <- plyr::mapvalues(seurat$SexChrom,
                                            from = c("XX", "XXY", "XY", "XYY"),
                                            to = c(2, 2, 1 ,1)) %>% as.character() %>% as.numeric()
seurat$Ydose <- plyr::mapvalues(seurat$SexChrom,
                                       from = c("XX", "XXY", "XY", "XYY"),
                                       to = c(0, 1, 1 ,2)) %>% factor(levels = c(0, 1, 2))
seurat$Ydose_cont <- plyr::mapvalues(seurat$SexChrom,
                                            from = c("XX", "XXY", "XY", "XYY"),
                                            to = c(0, 1, 1 ,2)) %>% as.character() %>% as.numeric()
seurat$base_SexChrom <- plyr::mapvalues(seurat$SexChrom,
                                               from = c("XX", "XXY", "XY", "XYY"),
                                               to = c("XX", "XY", "XY", "XY")) %>% factor(levels = c("XY", "XX"))
seurat$addX <- plyr::mapvalues(seurat$SexChrom,
                                      from = c("XX", "XXY", "XY", "XYY"),
                                      to = c(0, 1, 0, 0)) %>% as.numeric()
seurat$addY <- plyr::mapvalues(seurat$SexChrom,
                                      from = c("XX", "XXY", "XY", "XYY"),
                                      to = c(0, 0, 0, 1)) %>% as.numeric()

seurat$Gonad <- factor(seurat$Gonad, levels = c("Male", "Female"))
seurat$SexChrom <- factor(seurat$SexChrom, levels = c("XY", "XX", "XYY", "XXY"))
seurat$Genotype <- factor(seurat$Genotype, levels = c("XYM", "XYF", "XXM", "XXF", "XYYM", "XYYF", "XXYM", "XXYF"))
seurat$Trisomy <- factor(seurat$Trisomy, levels = c("Non-Trisomy", "Trisomy"))

############################################################
##Step 2 : Run Limma on both models 
############################################################

##normative analysis (XYM vs XXM) 
resultlist_normative <- list()
for (i in 1:length(allcells)) {
  currentcell <- as.character(allcells[i])
  currentsubset <- subset(seurat, cell.type == currentcell)
  
  dge0 <- DGEList(currentsubset@assays$RNA@counts, group = currentsubset$data.sampleName)
  dge <- dge0[!rowSums(dge0$counts==0)==ncol(dge0$counts), ]
  keep.exprs <- filterByExpr(dge)
  dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~Genotype, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3) 
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  resultlist_normative[["XXF_vs_XYM"]] <- topTable(fit, n = Inf, adjust.method = "BH", coef = "GenotypeXXF")
  save(resultlist_normative, file = paste0("~/scratch/medial_septum_sct/Data/Derived/Pseudobulk_DEG/normative/", currentcell,".rda"))
}

##sex effect analysis (Gonad, SexChrom, addX, addY) 
resultlist_sexeffects <- list()
for (i in 1:length(allcells)) {
  currentcell <- as.character(allcells[i])
  currentsubset <- subset(seurat, cell.type == currentcell)
  
  dge0 <- DGEList(currentsubset@assays$RNA@counts, group = currentsubset$data.sampleName)
  dge <- dge0[!rowSums(dge0$counts==0)==ncol(dge0$counts), ]
  keep.exprs <- filterByExpr(dge)
  dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~Gonad + base_SexChrom + addX + addY, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3) 
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  resultlist_sexeffects[["Ovary_vs_Testis"]] <- topTable(fit, n = Inf, adjust.method = "BH", coef = "GonadFemale")
  resultlist_sexeffects[["XX_vs_XY"]] <- topTable(fit, n = Inf, adjust.method = "BH", coef = "base_SexChromXX")
  resultlist_sexeffects[["addX"]] <- topTable(fit, n = Inf, adjust.method = "BH", coef = "addX1") ##check coefficient name - might be addX 
  resultlist_sexeffects[["addY"]] <- topTable(fit, n = Inf, adjust.method = "BH", coef = "addY1")
  
  save(resultlist_normative, file = paste0("~/scratch/medial_septum_sct/Data/Derived/Pseudobulk_DEG/sexeffects/", currentcell,".rda"))
}

############################################################
##Step 3 : Process DEGs into single dataframe 
############################################################

togenes <- list()
DEG_results <- list()
for(i in 1:length(allcells)) {
  
  currentcell <- as.character(allcells[i])
  
  load(paste0("~/scratch/medial_septum_sct/Data/Derived/Pseudobulk_DEG/normative/", currentcell,".rda"))
  load(paste0("~/scratch/medial_septum_sct/Data/Derived/Pseudobulk_DEG/sexeffects/", currentcell,".rda"))
  
  ##Normative Effect 
  for(k in 1:length(resultlist_normative)) {
    currentframe <- resultlist_normative[[k]]
    if(nrow(currentframe) == 0) {next}
    
    currentframe <- tibble::rownames_to_column(currentframe, "gene")
    currentcomparison <- names(resultlist_normative[k])
    currentframe$cell_type <- currentcell
    currentframe$comparison <- currentcomparison
    currentframe <- currentframe %>% dplyr::arrange(-logFC)
    
    topgenes[[currentcomparison]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
    DEG_results[[currentcomparison]] <- currentframe
  }
  
  ##Sex Effects 
  for(k in 1:length(resultlist_sexeffects)) {
    currentframe <- resultlist_sexeffects[[k]]
    if(nrow(currentframe) == 0) {next}
    
    currentframe <- tibble::rownames_to_column(currentframe, "gene")
    currentcomparison <- names(resultlist_sexeffects[k])
    currentframe$cell_type <- currentcell
    currentframe$comparison <- currentcomparison
    currentframe <- currentframe %>% dplyr::arrange(-logFC)
    
    topgenes[[currentcomparison]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
    DEG_results[[currentcomparison]] <- currentframe
  }
}

deg_df <- do.call("rbind", DEG_results)
topgenes_df <- do.call("rbind", topgenes)

save(deg_df, topgenes_df, file = "yourfilepath.rda")
