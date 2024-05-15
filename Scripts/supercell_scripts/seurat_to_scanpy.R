library(Seurat)
library(SeuratDisk)
library(magrittr)

setwd("~/scratch/medial_septum_sct")

seuratObject <- readRDS("Saves/For_Figures.rds")
seuratObject@active.assay <- "RNA"
seuratObject@meta.data[,'cell.type'] <- as.character(seuratObject@meta.data[,'cell.type'])
seuratObject@meta.data <- seuratObject@meta.data[,c('data.sampleName', 'Genotype', 'SexChrom', 'Gonad', 'cell.type', 'Trisomy')]
unlink("medial_septum.h5Seurat")
unlink("medial_septum.h5ad")
SaveH5Seurat(seuratObject, filename = "Saves/AnnData/medial_septum.h5Seurat")
Convert("Saves/AnnData/medial_septum.h5Seurat", dest = "h5ad")

#### #### #### #### #### #### #### 

##celltype loop - not needed
seuratObject <- readRDS("Saves/For_Figures.rds")
seuratObject@active.assay <- "RNA"
seurat_list <- SplitObject(seuratObject, split.by = "cell.type")

for (i in 1:length(seurat_list)) {
  
  seurat <- seurat_list[[i]]
  cell.type <- names(seurat_list)[i]
  data.dir <- paste0("Saves/AnnData/celltype_anndata/", cell.type, ".h5Seurerat")
  
  seurat@active.assay <- "RNA"
  seurat@meta.data[,'cell.type'] <- as.character(seurat@meta.data[,'cell.type'])
  seurat@meta.data <- seurat@meta.data[, c('data.sampleName', 'Genotype', 'SexChrom', 'Gonad', 'cell.type', 'Trisomy')]
  unlink(paste0(cell.type, ".h5Seurat"))
  unlink(paste0(cell.type, ".h5ad"))
  SaveH5Seurat(seurat, filename = data.dir)
  Convert(data.dir, dest = "h5ad")
  
}
