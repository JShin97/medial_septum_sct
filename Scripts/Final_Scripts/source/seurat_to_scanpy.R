library(Seurat)
library(SeuratDisk)
library(magrittr)

setwd("/u/scratch/j/jshin/medial_septum_sct")

##convert RNA assay 
seuratObject <- readRDS("Data/Derived/PostFilter/MS_cellbender_complete_analysis/seurat_harmony_sampleid_cellbender_df.Rds")
seuratObject@active.assay <- "RNA"
#seuratObject@meta.data[,'cell.type'] <- as.character(seuratObject@meta.data[,'cell.type'])
seuratObject@meta.data <- seuratObject@meta.data[,c('SampleID', 'Sample2', 'Sample_GT', 'Genotype', 'SexChrom', 'Gonad', 'Trisomy')]
unlink("MS_harmony_sampleid_Anndata.h5Seurat")
unlink("MS_harmony_sampleid_Anndata.h5ad")
SaveH5Seurat(seuratObject, filename = "MS_harmony_sampleid_Anndata.h5Seurat")
Convert("MS_harmony_sampleid_Anndata.h5Seurat", dest = "h5ad")

##convert SCT assay 
sct_counts <- seuratObject@assays$SCT@counts
sct_data <- seuratObject@assays$SCT@data

seuratObject2 <- CreateSeuratObject(counts = sct_counts, assay = "RNA")
seuratObject2@assays$RNA@data <- sct_data
if(identical(rownames(seuratObject2@meta.data), rownames(seuratObject@meta.data))) {seuratObject2@meta.data <- seuratObject@meta.data}
seuratObject2@meta.data <- seuratObject2@meta.data[,c('SampleID', 'Sample2', 'Sample_GT', 'Genotype', 'SexChrom', 'Gonad', 'Trisomy')]
unlink("MS_harmony_sampleid_SCT_Anndata.h5Seurat")
unlink("MS_harmony_sampleid_SCT_Anndata.h5ad")
SaveH5Seurat(seuratObject2, filename = "MS_harmony_sampleid_SCT_Anndata.h5Seurat")
Convert("MS_harmony_sampleid_SCT_Anndata.h5Seurat", dest = "h5ad")


#### #### #### #### #### #### #### 

##celltype loop 
seuratObject <- readRDS("/u/scratch/j/jshin/medial_septum_sct/Saves/MS_cellbender_complete_saves/seurat_harmony_sampleid_filtered_annotated.Rds")
seuratObject@active.assay <- "RNA"
seurat_list <- SplitObject(seuratObject, split.by = "celltype")

lapply(seq_along(seurat_list), function(i) {
  
  seurat <- seurat_list[[i]]
  cell.type <- names(seurat_list)[i]
  data.dir <- paste0("/u/scratch/j/jshin/medial_septum_sct/Data/Derived/network/SCING_complete/celltype_adata/", cell.type, ".h5Seurat")
  
  seurat@active.assay <- "RNA"
  seurat@meta.data[,'celltype'] <- as.character(seurat@meta.data[,'celltype'])
  seurat@meta.data <- seurat@meta.data[, c('SampleID', 'Genotype', 'SexChrom', 'Gonad', 'Trisomy', 'class_simplify', 'subclass_simplify', "neuronal", "NT", 'celltype')]
  unlink(paste0(cell.type, ".h5Seurat"))
  unlink(paste0(cell.type, ".h5ad"))
  SaveH5Seurat(seurat, filename = data.dir)
  Convert(data.dir, dest = "h5ad")
  
})

#### #### #### #### #### #### #### 

##celltype loop -- RNA counts only 
seuratObject <- readRDS("/u/scratch/j/jshin/medial_septum_sct/Saves/MS_cellbender_complete_saves/seurat_harmony_sampleid_filtered_annotated.Rds")
seuratObject@active.assay <- "RNA"
seurat_list <- SplitObject(seuratObject, split.by = "celltype")

lapply(seq_along(seurat_list), function(i) {
  
  seurat <- seurat_list[[i]]
  cell.type <- names(seurat_list)[i]
  data.dir <- paste0("/u/scratch/j/jshin/medial_septum_sct/Data/Derived/network/SCING_complete/celltype_RNA_adata/", cell.type, ".h5Seurat")
  
  counts <- seurat@assays$RNA@counts
  data <- seurat@assays$RNA@data
  
  seurat@meta.data[,'celltype'] <- as.character(seurat@meta.data[,'celltype'])
  metadata <- seurat@meta.data[, c('SampleID', 'Genotype', 'SexChrom', 'Gonad', 'Trisomy', 'class_simplify', 'subclass_simplify', "neuronal", "NT", 'celltype')]
  
  seurat_RNA <- CreateSeuratObject(counts = counts,
                                   meta.data = metadata,
                                   assay = "RNA")
  
  unlink(paste0(cell.type, ".h5Seurat"))
  unlink(paste0(cell.type, ".h5ad"))
  SaveH5Seurat(seurat_RNA, filename = data.dir)
  Convert(data.dir, dest = "h5ad")
  
})
