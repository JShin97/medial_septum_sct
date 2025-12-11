##DEG based on metacell counts

library(limma)
library(edgeR)
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(SingleCellExperiment)
library(ggpubr)
library(dplyr)
library(SuperCell)
library(hdWGCNA)

####################################################################
## Set metadata 
####################################################################
seurat <- readRDS("/u/home/j/jshin/scratch/medial_septum_sct/Saves/MS_cellbender_complete_saves/seurat_harmony_sampleid_filtered_annotated.Rds") %>% DietSeurat(counts = TRUE)

seurat$GenotypeaddX <- plyr::mapvalues(seurat$Genotype,
                                       from = c("XYM", "XYF", "XXYM", "XXYF"),
                                       to = c("XYM", "XYF", "XYM", "XYF")) %>% factor(levels = c("XYM", "XYF"))
seurat$addX <- plyr::mapvalues(seurat$SexChrom,
                               from = c("XY", "XXY"),
                               to = c(0, 1)) %>% as.numeric()
seurat$GenotypeaddY <- plyr::mapvalues(seurat$Genotype,
                                       from = c("XXM", "XXF", "XXYM", "XXYF", "XYM", "XYF", "XYYM", "XYYF"),
                                       to = c("XXM", "XXF", "XXM", "XXF", "XYM", "XYF", "XYM", "XYF")) %>% factor(levels = c("XYM", "XYF", "XXM", "XXF"))
seurat$addY <- plyr::mapvalues(seurat$SexChrom,
                               from = c("XX", "XXY", "XY", "XYY"),
                               to = c(0, 1, 0, 2)) %>% factor(., levels = c(0, 1, 2))
seurat$addY1 <- plyr::mapvalues(seurat$SexChrom,
                                from = c("XX", "XXY", "XY", "XYY"),
                                to = c(0, 1, 0, 0)) %>% as.character %>% as.numeric()
seurat$addY2 <- plyr::mapvalues(seurat$SexChrom,
                                from = c("XX", "XXY", "XY", "XYY"),
                                to = c(0, 0, 0, 1)) %>% as.character %>% as.numeric()

seurat$Gonad <- factor(seurat$Gonad, levels = c("Male", "Female"))
seurat$SexChrom <- factor(seurat$SexChrom, levels = c("XY", "XX", "XYY", "XXY"))
seurat$Genotype <- factor(seurat$Genotype, levels = c("XYM", "XYF", "XXM", "XXF", "XYYM", "XYYF", "XXYM", "XXYF"))
seurat$Trisomy <- factor(seurat$Trisomy, levels = c("Non-Trisomy", "Trisomy"))
seurat$SampleID <- gsub("-", "", seurat$SampleID)
seurat$Sequencing_Date <- factor(seurat$Sequencing_Date) %>% gsub("-", "_", .)

sample_list <- SplitObject(seurat, split.by = "SampleID")

####################################################################
## Metacell construction
####################################################################
gamma <- 20 # Graining level

##sample-wise construction 
makeSeuratMC <- function(seurat_object, 
                         genes = NULL, 
                         metaFields = c("SampleID", "celltype", "Genotype", "Gonad", "SexChrom", "Trisomy", "GenotypeaddX", "GenotypeaddY", "addX", "addY1", "addY2", "Sequencing_Date"), 
                         returnMC = T) {
  
  if(is.null(genes)) {
    genes <- VariableFeatures(seurat_object)
  }
  
  MC <- SCimplify(GetAssayData(seurat_object, slot = "data"), 
                  n.pc = 20, 
                  k.knn = 5, 
                  gamma = 20, 
                  genes.use = genes)
  
  MC$purity <- supercell_purity(clusters = seurat_object$celltype, 
                                supercell_membership = MC$membership) 
  
  for(m in metaFields) {
    MC[[m]] <- supercell_assign(clusters = seurat_object@meta.data[, m], 
                                supercell_membership = MC$membership, 
                                method = "absolute") 
  }
  
  GE <- supercell_GE(as.matrix(GetAssayData(seurat_object, slot = "data")),
                     groups = MC$membership)
  
  seuratMC <- supercell_2_Seurat(SC.GE = GE, MC, fields = c(metaFields, "purity"))
  
  res <- seuratMC
  if(returnMC) {
    res <- list(seuratMC = seuratMC, 
                SC = MC) 
  }
  return(res)
}

MC.list <- list()
for(i in 1:length(sample_list)) {
  st <- sample_list[[i]]
  sample <- unique(st$SampleID)
  
  seuratMC <- makeSeuratMC(st, returnMC = F)
  seuratMC <- subset(seuratMC, purity >= 0.7)
  
  MC.list[[sample]] <- seuratMC
}

features <- Seurat::SelectIntegrationFeatures(object.list = MC.list, nfeatures = 2000)
nSingleCells <- 0
for (i in 1:length(MC.list)) {
  MC.list[[i]] <- RenameCells(MC.list[[i]],add.cell.id = unique(MC.list[[i]]$SampleID))
  MC.list[[i]] <- ScaleData(MC.list[[i]],features = features)
  MC.list[[i]] <- RunPCA(MC.list[[i]] ,features = features,npcs = 20)
  nSingleCells <- nSingleCells + sum(MC.list[[i]]$size)
}
nSingleCells

integrated <- scCustomize::Merge_Seurat_List(MC.list)
DefaultAssay(integrated) <- "RNA"
integrated <- ScaleData(integrated, features = features)
integrated <- RunPCA(integrated, npcs = 20, verbose = FALSE,features = features)
integrated <- RunHarmony(integrated, c("Sequencing_Date"), theta = 1)
integrated <- RunUMAP(integrated, reduction = "harmony",
                      reduction.name = "umap.harmony",
                      dims = 1:20)

DimPlot(integrated, group.by = "SampleID", reduction = "umap.harmony")
DimPlot(integrated, group.by = "Genotype", reduction = "umap.harmony")
DimPlot(integrated, group.by = "Sequencing_Date", reduction = "umap.harmony")
DimPlot(integrated, group.by = "celltype", reduction = "umap.harmony")

Idents(integrated) <- "celltype"
scCustomize::Cluster_Highlight_Plot(integrated, cluster_name = "GABA-Chol", reduction = "umap.harmony")

saveRDS(integrated, "/u/home/j/jshin/scratch/medial_septum_sct/Saves/MS_cellbender_complete_saves/supercell_seurat_harmony_sampleid_filtered_annotated.Rds")

####################################################################
## DEG analysis with metacells
####################################################################
sc_seurat <- readRDS("/u/home/j/jshin/scratch/medial_septum_sct/Saves/MS_cellbender_complete_saves/supercell_seurat_harmony_sampleid_filtered_annotated.Rds")

sc_seurat_list <- SplitObject(sc_seurat, split.by = "celltype")
sc_seurat_list <- sc_seurat_list[!names(sc_seurat_list) %in% "Misc_NT"]
assay_to_use <- "RNA" 

##normative analysis (XYM vs XXM) 
for(i in seq_along(sc_seurat_list)) {
  
  resultlist_normative <- list()
  currentcell <- names(sc_seurat_list[i])
  currentsubset <- sc_seurat_list[[i]]
  
  gene_counts <- currentsubset[[assay_to_use]]@counts
  gene_data <- currentsubset[[assay_to_use]]@data
  
  genes_to_keep <- which(Matrix::rowMeans(gene_counts > 0) >= 0.15)
  gene_data <- gene_data[genes_to_keep,]
  
  design <- model.matrix(~0 + Genotype + Sequencing_Date, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- gene_data 
  fit <- lmFit(y, design = design)
  contrasts <- makeContrasts(XXF_vs_XYM = GenotypeXXF - GenotypeXYM, levels = design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  results <- topTable(fit, n = Inf, coef = "XXF_vs_XYM", adjust.method = "BH")
  resultlist_normative[["XXF_vs_XYM"]] <- results
  save(resultlist_normative, file = paste0("Results/complete_analysis/DEG_SD_metacell/normative/", currentcell,".rda"))
  
}

##sexchrom and gonad analysis (FCG)

for(i in seq_along(sc_seurat_list)) {
  
  resultlist_fcg <- list()
  currentcell <- names(sc_seurat_list[i])
  currentsubset <- sc_seurat_list[[i]]
  currentsubset <- subset(currentsubset, Genotype %in% c("XXM", "XYM", "XXF", "XYF"))
  currentsubset$Genotype <- factor(currentsubset$Genotype, levels = c("XYM", "XYF", "XXM", "XXF"))
  currentsubset$SexChrom <- factor(currentsubset$SexChrom, levels = c("XY", "XX"))
  currentsubset$Gonad <- factor(currentsubset$Gonad, levels = c("Male", "Female"))
  
  gene_counts <- currentsubset[[assay_to_use]]@counts
  gene_data <- currentsubset[[assay_to_use]]@data
  
  genes_to_keep <- which(Matrix::rowMeans(gene_counts > 0) >= 0.15) 
  gene_data <- gene_data[genes_to_keep,]
  
  design <- model.matrix(~Gonad + SexChrom + Gonad:SexChrom + Sequencing_Date, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- gene_data 
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  resultlist_fcg[["Ovary_vs_Testis"]] <- topTable(fit, n = Inf, coef = "GonadFemale", adjust.method = "BH")
  resultlist_fcg[["XX_vs_XY"]] <- topTable(fit, n = Inf, coef = "SexChromXX", adjust.method = "BH")
  resultlist_fcg[["Ovary:SexChromXX"]] <- topTable(fit, n = Inf, coef = "GonadFemale:SexChromXX", adjust.method = "BH")
  
  save(resultlist_fcg, file = paste0("Results/complete_analysis/DEG_SD_metacell/fcg/", currentcell,".rda"))
}

##addX analysis 
for(i in seq_along(sc_seurat_list)) {
  
  resultlist_addx <- list()
  currentcell <- names(sc_seurat_list[i])
  currentsubset <- sc_seurat_list[[i]]
  currentsubset <- subset(currentsubset, Genotype %in% c("XYM", "XXYM", "XYF", "XXYF"))
  currentsubset$addX <- as.numeric(currentsubset$addX)
  
  gene_counts <- currentsubset[[assay_to_use]]@counts
  gene_data <- currentsubset[[assay_to_use]]@data
  
  genes_to_keep <- which(Matrix::rowMeans(gene_counts > 0) >= 0.15)
  gene_data <- gene_data[genes_to_keep,]
  
  design <- model.matrix(~GenotypeaddX + addX + Sequencing_Date, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- gene_data 
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  resultlist_addx[["addX"]] <- topTable(fit, n = Inf, coef = "addX", adjust.method = "BH")
  
  save(resultlist_addx, file = paste0("Results/complete_analysis/DEG_SD_metacell/addx/", currentcell,".rda"))
}

##addY1 analysis 
for(i in seq_along(sc_seurat_list)) {
  
  resultlist_addy1 <- list()
  currentcell <- names(sc_seurat_list[i])
  currentsubset <- sc_seurat_list[[i]]
  currentsubset <- subset(currentsubset, Genotype %in% c("XXM", "XXYM", "XXF", "XXYF"))
  currentsubset$addY1 <- as.numeric(currentsubset$addY1)
  
  gene_counts <- currentsubset[[assay_to_use]]@counts
  gene_data <- currentsubset[[assay_to_use]]@data
  
  genes_to_keep <- which(Matrix::rowMeans(gene_counts > 0) >= 0.15)  
  gene_data <- gene_data[genes_to_keep,]
  
  design <- model.matrix(~GenotypeaddY + addY1 + Sequencing_Date, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- gene_data 
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  resultlist_addy1[["addY1"]] <- topTable(fit, n = Inf, coef = "addY1", adjust.method = "BH")
  
  save(resultlist_addy1, file = paste0("Results/complete_analysis/DEG_SD_metacell/addy1/", currentcell,".rda"))
}

##addY2 analysis 
for(i in seq_along(sc_seurat_list)) {
  
  resultlist_addy2 <- list()
  currentcell <- names(sc_seurat_list[i])
  currentsubset <- sc_seurat_list[[i]]
  currentsubset <- subset(currentsubset, Genotype %in% c("XYM", "XYYM", "XYF", "XYYF"))
  currentsubset$addY2 <- as.numeric(currentsubset$addY2)
  
  gene_counts <- currentsubset[[assay_to_use]]@counts
  gene_data <- currentsubset[[assay_to_use]]@data
  
  genes_to_keep <- which(Matrix::rowMeans(gene_counts > 0) >= 0.15)
  gene_data <- gene_data[genes_to_keep,]
  
  design <- model.matrix(~GenotypeaddY + addY2 + Sequencing_Date, data = currentsubset@meta.data)
  
  y <- new("EList")
  y$E <- gene_data 
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  resultlist_addy2[["addY2"]] <- topTable(fit, n = Inf, coef = "addY2", adjust.method = "BH")
  
  save(resultlist_addy2, file = paste0("Results/complete_analysis/DEG_SD_metacell/addy2/", currentcell,".rda"))
}

####################################################################
## DEG Processing
####################################################################

data.dir <- "Results/complete_analysis/DEG_SD_metacell/"
allfiles <- list.files(data.dir)
allcells <- list.files(paste0(data.dir, allfiles[1])) %>% gsub(".rda", "", .)

topgenes <- list()
DEG_results <- list()

for(i in 1:length(allcells)) {
  currentcell <- allcells[i]
  
  #Normative difference 
  file.dir <- paste0(data.dir, "normative/", currentcell, ".rda")
  load(file.dir)
  
  for(k in 1:length(resultlist_normative)) {
    
    currentframe <- resultlist_normative[[k]]
    if(nrow(currentframe) == 0) {next}
    
    currentframe <- tibble::rownames_to_column(currentframe, "gene")
    currentcomparison <- names(resultlist_normative[k])
    currentframe$cell_type <- currentcell
    currentframe$comparison <- currentcomparison
    currentframe <- currentframe %>% dplyr::arrange(-logFC)
    currentmodel <- "Normative"
    
    topgenes[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
    DEG_results[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe
  }
  
  #FCG
  file.dir <- paste0(data.dir, "fcg/", currentcell, ".rda")
  load(file.dir)
  
  for(k in 1:length(resultlist_fcg)) {
    
    currentframe <- resultlist_fcg[[k]]
    if(nrow(currentframe) == 0) {next}
    
    currentframe <- tibble::rownames_to_column(currentframe, "gene")
    currentcomparison <- names(resultlist_fcg[k])
    currentframe$cell_type <- currentcell
    currentframe$comparison <- currentcomparison
    currentframe <- currentframe %>% dplyr::arrange(-logFC)
    currentmodel <- "FCG"
    
    topgenes[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
    DEG_results[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe
  }
  
  #addX
  file.dir <- paste0(data.dir, "addx/", currentcell, ".rda")
  load(file.dir)
  
  currentframe <- resultlist_addx[["addX"]]
  if(nrow(currentframe) == 0) {next}
  
  currentframe <- tibble::rownames_to_column(currentframe, "gene")
  currentcomparison <- names(resultlist_addx["addX"])
  currentframe$cell_type <- currentcell
  currentframe$comparison <- currentcomparison
  currentframe <- currentframe %>% dplyr::arrange(-logFC)
  currentmodel <- "addX"
  
  topgenes[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
  DEG_results[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe
  
  #addY1
  file.dir <- paste0(data.dir, "addy1/", currentcell, ".rda")
  load(file.dir)
  
  currentframe <- resultlist_addy1[["addY1"]]
  if(nrow(currentframe) == 0) {next}
  
  currentframe <- tibble::rownames_to_column(currentframe, "gene")
  currentcomparison <- names(resultlist_addy1["addY1"])
  currentframe$cell_type <- currentcell
  currentframe$comparison <- currentcomparison
  currentframe <- currentframe %>% dplyr::arrange(-logFC)
  currentmodel <- "addY1"
  
  topgenes[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
  DEG_results[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe
  
  #addY2
  file.dir <- paste0(data.dir, "addy2/", currentcell, ".rda")
  load(file.dir)
  
  currentframe <- resultlist_addy2[["addY2"]]
  if(nrow(currentframe) == 0) {next}
  
  currentframe <- tibble::rownames_to_column(currentframe, "gene")
  currentcomparison <- names(resultlist_addy2["addY2"])
  currentframe$cell_type <- currentcell
  currentframe$comparison <- currentcomparison
  currentframe <- currentframe %>% dplyr::arrange(-logFC)
  currentmodel <- "addY2"
  
  topgenes[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe$gene[c(1:min(3,nrow(currentframe)))]
  DEG_results[[currentmodel]][[currentcomparison]][[currentcell]] <- currentframe
}

DEG_df <- do.call("rbind", unlist(DEG_results, recursive=FALSE)) %>% bind_rows()
DEG_df$effect <- plyr::mapvalues(DEG_df$comparison, 
                                 from = c("XXF_vs_XYM", "Ovary_vs_Testis", "XX_vs_XY", "Ovary:SexChromXX", "addX", "addY1", "addY2"),
                                 to = c("Typical", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2"))

save(topgenes, DEG_results, file = "Results/complete_analysis/DEG_processed/DEG_SD_metacell/topgenes_DEGlist_RNA_SD.rda")
save(DEG_df, file = "Results/complete_analysis/DEG_processed/DEG_SD_metacell/DEG_df_RNA_SD.rda")
