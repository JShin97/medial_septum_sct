##Script for all QC Scripts 

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = "/u/project/xyang123/jshin/medial_septum_sct/")

library(Seurat)
library(stringr)
library(readxl)
library(tidyr)
library(plyr)
library(dplyr)
library(readr)
library(multtest)
library(mutoss)
library(metap)
library(purrr)
library(SoupX)
library(hdf5r)
library(ggpubr)
library(RSpectra)
library(DropletUtils)
library(BiocParallel)
library(DoubletFinder)
library(harmony)

source("./source/scFunctions.R")

####################################################################
##Part 1 : Load cellbender processed raw data 
####################################################################

sample.dir <- "/u/home/j/jshin/scratch/medial_septum_sct/Data/Raw/Cellbender_counts2"
samples <- list.dirs(path = sample.dir, full.names = FALSE, recursive = FALSE)

##load original objects
seurat_list <- lapply(samples, function(sample){
  data_path <- paste0(sample.dir, "/", sample, "/", sample, "_filtered_seurat.h5")
  count_data <- Read10X_h5(data_path)
  cur_seurat <- CreateSeuratObject(
    counts = count_data,
    project='MS_SCT'
  )
  cur_seurat$SampleID <- sample
  return(cur_seurat)
})
saveRDS(seurat_list, "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PreFilter/seurat_list_cellbender.Rds")

####################################################################
##Part 2 : Run DoubletFinder
####################################################################
seurat_list <- readRDS("/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PreFilter/seurat_list_cellbender.Rds")

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
  png(file = paste0("Plots/QC/doubletfinder/", unique(sample$SampleID), ".png"),width=3,height=5,units="in",res=600)
  bcmvn_sample <- find.pK(sweep.stats)
  dev.off()
  mpK <-as.numeric(as.vector(bcmvn_sample$pK[which.max(bcmvn_sample$BCmetric)]))
  # (3). Run DoubletFinder
  DoubletRate = ncol(sample)*8*1e-6 # ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
  nExp_poi <- round(DoubletRate*ncol(sample))
  sample <- doubletFinder(sample, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  # (4). Add doubletfinder metadata
  table(sample@meta.data[,paste0("DF.classifications_0.25_",mpK,"_",nExp_poi)])
  sample$doubletfinder<-sample@meta.data[,paste0("DF.classifications_0.25_",mpK,"_",nExp_poi)]
  
  return(sample)
})

saveRDS(seurat_list_df, file = "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PreFilter/seurat_list_cellbender_df.Rds")

####################################################################
##Part 3 : Add metadata 
####################################################################
seuratObject <- merge(x = seurat_list_df[[1]], y = seurat_list_df[2:length(seurat_list_df)], merge.data=TRUE, add.cell.ids = samples)

seuratObject@meta.data[, grepl("pANN", colnames(seuratObject@meta.data))] <- NULL
seuratObject@meta.data[, grepl("DF", colnames(seuratObject@meta.data))] <- NULL

##add metadata
load(file = "Data/metadata/metadata.Rda")
sample_metadata <- seuratObject@meta.data
sample_metadata <- sample_metadata %>%
  tibble::rownames_to_column(var = "cell") %>%
  left_join(., original_metadata, by = c("SampleID" = "orig.ident")) %>%
  left_join(., technical_metadata, by = c("SampleID" = "Sample_ID"))
sample_metadata <- sample_metadata %>%
  dplyr::select(-Genotype.y) %>%
  dplyr::rename(Genotype = Genotype.x)
rownames(sample_metadata) <- sample_metadata$cell
sample_metadata <- sample_metadata[, -1]
for (meta in names(sample_metadata)) {
  seuratObject[[meta]] <- sample_metadata[[meta]]
}

##remap metadata based on validation results
seuratObject$Genotype_original <- seuratObject$Genotype

metadata <- seuratObject@meta.data
metadata[metadata$SampleID == "TriSeptum-2", "Genotype"] <- "XXYF"
metadata[metadata$SampleID == "TriSeptum-5", "Genotype"] <- "XXF"
metadata[metadata$SampleID == "TriSeptum-12", "Genotype"] <- "XYM"
metadata[metadata$SampleID == "TriSeptum-24", "Genotype"] <- "XYM"

if(identical(row.names(seuratObject@meta.data), row.names(metadata))) {seuratObject@meta.data <- metadata}

saveRDS(seuratObject, file = "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PreFilter/seurat_merged_cellbender_df.Rds")

####################################################################
##Part 4 : Filtering 
####################################################################
seuratObject = getCellQuality(seuratObject = seuratObject,
                              feature_patterns=list("percent.mito"="^mt-",
                                                    "percent.ribo"=c("^Rps","^Rpl"),
                                                    "percent.pred"=c("^Gm1","^Gm2","^Gm3","^Gm4","^Gm5","^Gm6","^Gm7","^Gm8","^Gm9"),
                                                    "percent.Hb"=c("^Hba-","^Hbb-","^Hbq")))
cellQualityPlot(seuratObject=seuratObject,
                fileName="Plots/QC/ViolinPlots/MS_cellbender_PreFilter.pdf",
                featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito",
                               "percent.ribo","percent.pred","percent.Hb"),
                identPlot = "SampleID",
                H=9,W=20,
                pointSize=0)

filtered_sample = subset(x = seuratObject, 
                         subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 7500 & 
                           nCount_RNA > 500 &
                           nCount_RNA < 35000 &
                           percent.mito < .03 & 
                           percent.pred < .10 & 
                           percent.Hb < 0.05)

# redo cell quality plot to see QC after filtering
cellQualityPlot(seuratObject=filtered_sample,
                fileName=paste0("Plots/QC/ViolinPlots/MS_cellbender_PostFilter.pdf"),
                H=9,W=20,
                featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito",
                               "percent.ribo","percent.pred","percent.Hb"),
                identPlot = "SampleID",
                pointSize=0)

####################################################################
##Part 5 : Preliminary Clustering 
####################################################################
filtered_sample <- NormalizeData(filtered_sample, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_sample <- FindVariableFeatures(filtered_sample, selection.method = "vst", nfeatures = 2000)
filtered_sample <- ScaleData(filtered_sample, features = rownames(filtered_sample))
filtered_sample <- RunPCA(filtered_sample, features = VariableFeatures(object = filtered_sample))
filtered_sample <- FindNeighbors(filtered_sample, dims = 1:10)
filtered_sample <- RunUMAP(filtered_sample, dims = 1:10)
filtered_sample <- FindClusters(filtered_sample, resolution = c(0.5, 0.8, 1.0))

saveRDS(filtered_sample, file = "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_merged_cellbender_df.Rds")

##visualization 
load(file = "Data/metadata/metadata.Rda")
seuratObject <- readRDS("/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_merged_cellbender_df.Rds")

metadata <- seuratObject@meta.data
technical_metadata <- technical_metadata %>% group_by(Genotype) %>% mutate(id = row_number()) %>% ungroup() %>% dplyr::select(Sample_ID, id)
metadata <- metadata %>%
  rownames_to_column(var = "cellID") %>%
  left_join(., technical_metadata, by = c("SampleID" = "Sample_ID"))
rownames(metadata) <- metadata$cellID 
metadata$cellID <- NULL

##Add Sample2 and Sample_GT 
metadata$Sample2 <- paste0(metadata$Genotype, metadata$id)
metadata$Sample_GT <- paste(metadata$SampleID, metadata$Genotype, sep = "_")
metadata$Genotype <- factor(metadata$Genotype, levels = c("XYM", "XYF", "XXM", "XXF", "XYYM", "XYYF", "XXYM", "XXYF"))
sample2_order <- metadata$Sample2[order(metadata$Genotype, metadata$id)]
metadata$Sample2 <- factor(metadata$Sample2, levels = unique(sample2_order))
samplegt_order <- metadata$Sample_GT[order(metadata$Sample2)]
metadata$Sample_GT <- factor(metadata$Sample_GT, levels = unique(samplegt_order))

if(identical(row.names(seuratObject@meta.data), row.names(metadata))) {seuratObject@meta.data <- metadata}

seurat_sample_list <- SplitObject(seuratObject, split.by = "Sample_GT")
seurat_sample_list <- seurat_sample_list[levels(seuratObject$Sample_GT)]
umap_list <- lapply(seq_along(seurat_sample_list), function(i) {
  sample <- seurat_sample_list[[i]]
  plot <- DimPlot(sample, group.by = "seurat_clusters", label = TRUE) + ggtitle(as.character(names(seurat_sample_list[i])))
  return(plot)
})

png("Plots/QC/UMAP_preintegration/sample_umap.png", height = 30, width = 40, units = "in", res = 600)
ggarrange(plotlist = umap_list, ncol = 6, nrow = 4)
dev.off()

metadata <- seuratObject@meta.data
metadata$Sequencing_Date <- factor(seuratObject$Sequencing_Date, levels = unique(seuratObject$Sequencing_Date)[order(as.Date(unique(seuratObject$Sequencing_Date)))])
sample_sd_order <- metadata$Sample_GT[order(metadata$Sequencing_Date)]

seurat_sample_list <- seurat_sample_list[unique(sample_sd_order)]
umap_list <- lapply(seq_along(seurat_sample_list), function(i) {
  sample <- seurat_sample_list[[i]]
  sample_name <- as.character(names(seurat_sample_list[i]))
  plot <- DimPlot(sample, group.by = "seurat_clusters", label = TRUE) + ggtitle(paste(sample_name, unique(sample$Sequencing_Date), sep = "_"))
  return(plot)
})

png("Plots/QC/UMAP_preintegration/sample_sd_umap.png", height = 30, width = 40, units = "in", res = 600)
ggarrange(plotlist = umap_list, ncol = 6, nrow = 4)
dev.off()

####################################################################
##Part 6 : Integration  
####################################################################
seuratObject <- readRDS("/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_merged_cellbender_df.Rds")

##Integration by sample -- use these results 
seurat_sample_list <- SplitObject(seuratObject, split.by = "SampleID")

seurat_sample_list <- lapply(X = seurat_sample_list,
                             FUN = SCTransform,
                             vst.flavor = "v2",
                             method = "glmGamPoi",
                             return.only.var.genes = FALSE,
                             conserve.memory = TRUE)
var.features <- SelectIntegrationFeatures(object.list = seurat_sample_list, nfeatures = 3000)

seuratObject.sct <- merge(x = seurat_sample_list[[1]], y = seurat_sample_list[2:length(seurat_sample_list)], merge.data=TRUE)
VariableFeatures(seuratObject.sct) <- var.features
seuratObject.sct <- RunPCA(seuratObject.sct, verbose = TRUE)
seuratObject.sct <- RunHarmony(seuratObject.sct, assay.use="SCT", group.by.vars = "SampleID")
seuratObject.sct <- RunUMAP(seuratObject.sct, reduction = "harmony", dims = 1:30)
seuratObject.sct <- FindNeighbors(seuratObject.sct, reduction = "harmony", dims = 1:30)
seuratObject.sct <- FindClusters(seuratObject.sct, reduction = "harmony", resolution = c(0.5, 0.8, 1.0))

##set metadata
seuratObject.sct$Sample2 <- factor(seuratObject.sct$Sample2, levels = levels(seuratObject$Sample2))
seuratObject.sct$Sample_GT <- factor(seuratObject.sct$Sample_GT, levels = levels(seuratObject$Sample_GT))

metadata <- seuratObject.sct@meta.data
metadata$Sequencing_Date <- factor(seuratObject.sct$Sequencing_Date, levels = unique(seuratObject.sct$Sequencing_Date)[order(as.Date(unique(seuratObject.sct$Sequencing_Date)))])
sample_sd_order <- metadata$Sample_GT[order(metadata$Sequencing_Date)]

saveRDS(seuratObject.sct, file = "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_harmony_sampleid_cellbender_df.Rds")

####################################################################
##Part 7 : Prep for MapMyCells -- use seurat_to_scanpy.R script
##Part 7.5 : Run MapMyCells
####################################################################

####################################################################
##Part 8 : Process MapMyCells results
####################################################################
seuratObject.sct <- readRDS("/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_harmony_sampleid_cellbender_df.Rds")

##Load MapMyCells results
mapping <- read.csv(Sys.glob("Results/MapMyCells/MS_harmony_sampleid_SCT_Anndata_MapMyCells_10xWholeMouseBrain_HierarchicalMapping/*.csv"), comment.char="#")

data.frame(Cell_counts = sort(table(mapping$class_name), decreasing = T)) 
data.frame(Cell_counts = sort(table(mapping$subclass_name), decreasing = T)) 

cellcounts <- table(mapping$class_name) %>% as.data.frame()
subcellcounts <- table(mapping$subclass_name) %>% as.data.frame()

#Assign rare classes and subclasses as "other"
mapping$class_new <- mapping$class_name
mapping$class_new[mapping$class_name %in% cellcounts[cellcounts$Freq < 100, "Var1"]] <- "Other"
mapping$subclass_new <- mapping$subclass_name
mapping$subclass_new[mapping$subclass_name %in% subcellcounts[subcellcounts$Freq < 100, "Var1"]] <- "Other"

#Simplify class metadata -- neuronal vs non-neuronal
mapping$neuronal <- "Neuronal"
mapping[grepl("NN", mapping$subclass_new), "neuronal"] <- "Non-neuronal"

#Simplify by Neuron NT type 
mapping$NT <- mapping$neuronal
mapping[grepl("Gaba", mapping$subclass_new), "NT"] <- "GABA"
mapping[grepl("Glut", mapping$subclass_new), "NT"] <- "Glut"
#mapping[grepl("Gaba-Glut", mapping$subclass_new), "NT"] <- "GABA-Glut"
#mapping[grepl("Glyc", mapping$subclass_new), "NT"] <- "Glyc-GABA"
#mapping[grepl("Chol", mapping$subclass_new), "NT"] <- "Chol"
mapping[grepl("Gaba-Chol", mapping$subclass_new), "NT"] <- "GABA-Chol"
#mapping[grepl("Dopa", mapping$subclass_new), "NT"] <- "Dopa"
#mapping[grepl("Dopa-Gaba", mapping$subclass_new), "NT"] <- "Dopa-GABA"
#mapping[grepl("Sero", mapping$subclass_new), "NT"] <- "Glut-Sero"
#mapping[grepl("Hist", mapping$subclass_new), "NT"] <- "Hist-GABA"
mapping[grepl("Inh", mapping$subclass_new), "NT"] <- "Inh_IMN"
mapping[grepl("Neuronal", mapping$NT), "NT"] <- "Misc_NT" ##classify the rest as misc 

#Simplify by broad class 
mapping$class_simplify <- mapping$NT
mapping[mapping$NT == "Non-neuronal", "class_simplify"] <- mapping[mapping$NT == "Non-neuronal", "class_new"]
mapping$class_simplify <- plyr::mapvalues(mapping$class_simplify,
                                          from = c("30 Astro-Epen", "31 OPC-Oligo", "33 Vascular", "34 Immune"),
                                          to = c("Astro-Epen", "OPC-Oligo", "Vascular", "Immune")
)

#Simplify by subclass
mapping$subclass_simplify <- mapping$NT
mapping[mapping$NT == "Non-neuronal", "subclass_simplify"] <- mapping[mapping$NT == "Non-neuronal", "subclass_new"]
mapping$subclass_simplify <- plyr::mapvalues(mapping$subclass_simplify,
                                             from = c("318 Astro-NT NN", "319 Astro-TE NN", "323 Ependymal NN", 
                                                      "326 OPC NN", "327 Oligo NN", "330 VLMC NN", "331 Peri NN", "333 Endo NN", "334 Microglia NN"),
                                             to = c("Astro-NT", "Astro-TE", "Ependymal", "OPC", "Oligo", "VLMC", "Pericyte", 
                                                    "Endothelial", "Microglia"))

#Put row.names as data colnames and the order to match the data
rownames(mapping) <- mapping$cell_id
mapping <- mapping[rownames(seuratObject.sct@meta.data),]

#Add new metadata columns from mapping 
seuratObject.sct <- AddMetaData(seuratObject.sct,
                                metadata = mapping,
                                col.name = colnames(mapping))

####################################################################
##Part 9 : Annotate doublet clusters
####################################################################
seuratObject.sct$doubletcluster <- "FALSE"

##subcluster specific clusters 
seurat_subset <- subset(seuratObject.sct, SCT_snn_res.1 %in% c("20"))
DefaultAssay(seurat_subset) <- "RNA"
seurat_subset <- FindVariableFeatures(seurat_subset)
seurat_subset <- ScaleData(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- FindNeighbors(seurat_subset, dims = 1:50)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:50)
seurat_subset <- FindClusters(seurat_subset, resolution = c(0.03))
plot1 <- DimPlot(seurat_subset, group.by = c("doubletfinder", "seurat_clusters"))
plot2 <- DotPlot(seurat_subset, features = c("Hexb", "Tmem119", "C1qa", "Trem2", "Cx3cr1"))
ggsave(plot = plot1 + plot2, filename = "Plots/QC/doubletcluster/cluster20_subset.png", units = "in", height = 3, width = 12)
seurat_subset$doubletcluster <- ifelse(seurat_subset$RNA_snn_res.0.03 %in% c("1"), "TRUE", "FALSE")
seuratObject.sct@meta.data[rownames(seuratObject.sct@meta.data) %in% WhichCells(seurat_subset), "doubletcluster"] <- seurat_subset$doubletcluster

DimPlot(seuratObject.sct, group.by = c("doubletcluster", "class_simplify"))

seurat_subset <- subset(seuratObject.sct, SCT_snn_res.1 %in% c("9"))
DefaultAssay(seurat_subset) <- "RNA"
seurat_subset <- FindVariableFeatures(seurat_subset)
seurat_subset <- ScaleData(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- FindNeighbors(seurat_subset, dims = 1:50)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:50)
seurat_subset <- FindClusters(seurat_subset, resolution = c(0.05))
plot1 <- DimPlot(seurat_subset, group.by = c("doubletfinder", "seurat_clusters"))
ggsave(plot = plot1, filename = "Plots/QC/doubletcluster/cluster9_subset.png", units = "in", height = 3, width = 6)
seurat_subset$doubletcluster <- ifelse(seurat_subset$RNA_snn_res.0.05 %in% c("1", "2"), "TRUE", "FALSE")
seuratObject.sct@meta.data[rownames(seuratObject.sct@meta.data) %in% WhichCells(seurat_subset), "doubletcluster"] <- seurat_subset$doubletcluster

DimPlot(seuratObject.sct, group.by = c("doubletcluster", "class_simplify"))

seurat_subset <- subset(seuratObject.sct, SCT_snn_res.1 %in% c("18", "7", "13"))
DefaultAssay(seurat_subset) <- "RNA"
seurat_subset <- FindVariableFeatures(seurat_subset)
seurat_subset <- ScaleData(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- FindNeighbors(seurat_subset, dims = 1:50)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:50)
seurat_subset <- FindClusters(seurat_subset, resolution = c(0.05))
plot1 <- DimPlot(seurat_subset, group.by = c("doubletfinder", "seurat_clusters"))
ggsave(plot = plot1, filename = "Plots/QC/doubletcluster/cluster18_7_13_subset.png", units = "in", height = 3, width = 6)
seurat_subset$doubletcluster <- ifelse(seurat_subset$RNA_snn_res.0.05 %in% c("1", "2"), "TRUE", "FALSE")
seuratObject.sct@meta.data[rownames(seuratObject.sct@meta.data) %in% WhichCells(seurat_subset), "doubletcluster"] <- seurat_subset$doubletcluster

DimPlot(seuratObject.sct, group.by = c("doubletcluster", "class_simplify"))

seurat_subset <- subset(seuratObject.sct, SCT_snn_res.1 %in% c("2", "21", "38", "37", "35", "39"))
DefaultAssay(seurat_subset) <- "RNA"
seurat_subset <- FindVariableFeatures(seurat_subset)
seurat_subset <- ScaleData(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- RunPCA(seurat_subset, features = VariableFeatures(object = seurat_subset))
seurat_subset <- FindNeighbors(seurat_subset, dims = 1:50)
seurat_subset <- RunUMAP(seurat_subset, dims = 1:50)
seurat_subset <- FindClusters(seurat_subset, resolution = c(0.03))
plot1 <- DimPlot(seurat_subset, group.by = c("doubletfinder", "seurat_clusters"))
ggsave(plot = plot1, filename = "Plots/QC/doubletcluster/clusteroligo_subset.png", units = "in", height = 3, width = 6)
seurat_subset$doubletcluster <- ifelse(seurat_subset$RNA_snn_res.0.05 %in% c("1", "2"), "TRUE", "FALSE")
seuratObject.sct@meta.data[rownames(seuratObject.sct@meta.data) %in% WhichCells(seurat_subset), "doubletcluster"] <- seurat_subset$doubletcluster

saveRDS(seuratObject.sct, file = "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_harmony_sampleid_cellbender_doubletcluster.Rds")

####################################################################
##Part 10 : Create final saved object
####################################################################
seuratObject.sct <- readRDS(file = "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_harmony_sampleid_cellbender_doubletcluster.Rds")

seuratObject_filtered <- subset(seuratObject.sct, doubletcluster == "FALSE")

metadata <- seuratObject_filtered@meta.data
metadata$neuronal <- ifelse(metadata$SCT_snn_res.1 %in% c("2","21","37","35","39","20","33","9","36","25","34","27","18","7","13","36"),
                            "Non-neuronal","Neuronal")
metadata$celltype <- "NA"
metadata$celltype <- ifelse(metadata$SCT_snn_res.1 %in% c("0","28","10","29","24","16","19","8","1","6","22","15","4","14"), "GABA",
                            ifelse(metadata$SCT_snn_res.1 %in% c("23","30","5","26","17","31"), "Glut",
                                   ifelse(metadata$SCT_snn_res.1 %in% c("27"), "Inh_IMN",
                                          ifelse(metadata$SCT_snn_res.1 %in% c("2","21","35","39","33","37"), "Oligo",
                                                 ifelse(metadata$SCT_snn_res.1 %in% c("20"), "Microglia",
                                                        ifelse(metadata$SCT_snn_res.1 %in% c("9"), "OPC",
                                                               ifelse(metadata$SCT_snn_res.1 %in% c("25"), "Vascular", 
                                                                      ifelse(metadata$SCT_snn_res.1 %in% c("34"), "Ependymal",
                                                                             ifelse(metadata$SCT_snn_res.1 %in% c("18","7","13"), "Astrocyte",
                                                                                    metadata$subclass_simplify)))))))))
seuratObject_filtered$celltype <- metadata$celltype
seuratObject_filtered$neuronal <- metadata$neuronal
seuratObject_filtered$celltype <- ifelse(seuratObject_filtered$celltype %in% c("Astro-TE", "Astro-NT"), "Astrocyte", seuratObject_filtered$celltype)
seuratObject_filtered$celltype <- ifelse(seuratObject_filtered$celltype %in% c("Pericyte", "VLMC", "Endothelial"), "Vascular", seuratObject_filtered$celltype)

##fix metadata 
metadata <- seuratObject_filtered@meta.data
metadata$Gonad <- plyr::mapvalues(metadata$Genotype, 
                                  from = c("XYM", "XYF", "XXM", "XXF", "XYYM", "XYYF", "XXYM", "XXYF"), 
                                  to = c("Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female")) 
metadata$SexChrom <- plyr::mapvalues(metadata$Genotype, 
                                     from = c("XYM", "XYF", "XXM", "XXF", "XYYM", "XYYF", "XXYM", "XXYF"), 
                                     to = c("XY", "XY", "XX", "XX", "XYY", "XYY", "XXY", "XXY")) 
metadata$Trisomy <- plyr::mapvalues(metadata$Genotype, 
                                    from = c("XYM", "XYF", "XXM", "XXF", "XYYM", "XYYF", "XXYM", "XXYF"), 
                                    to = c("Non-Trisomy", "Non-Trisomy", "Non-Trisomy", "Non-Trisomy", "Trisomy", "Trisomy", "Trisomy", "Trisomy")) 

seuratObject_filtered@meta.data <- metadata 

saveRDS(seuratObject_filtered, file = "/u/home/j/jshin/scratch/medial_septum_sct/Saves/MS_cellbender_complete_saves/seurat_harmony_sampleid_filtered_annotated.Rds")