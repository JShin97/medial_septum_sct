library(Seurat)
library(SeuratDisk)
library(rhdf5)

setwd('~/scratch/medial_septum_sct')

ct <- c("Astrocyte", "Ependymal", "GABA-Chol", "GABA", "Glut", "Inh_IMN", "Microglia", "Misc_NT", "OPC", "Oligo", "Vascular")

for(i in 1:length(ct)) {
  
  adata_path <- paste0('Data/Derived/network/SCING_complete/supercells/', ct[i])
  h5ad_file <- paste0(adata_path, ".h5ad")
  
  seurat <- MuDataSeurat::ReadH5AD(h5ad_file)
  saveRDS(seurat, paste0('Data/Derived/network/SCING_complete/supercells_seurat/', ct[i], ".RDS"))
  
}

# for(i in 1:length(ct)) {
# 
#   adata_path <- paste0('Data/Derived/network/SCING_complete/supercells/', ct[i])
#   h5ad_file <- paste0(adata_path, ".h5ad")
#   Convert(h5ad_file, dest = "h5seurat", overwrite = T)
# 
#   seuratObj <- LoadH5Seurat(paste0(adata_path, ".h5seurat"), meta.data = FALSE, misc = FALSE)
# 
#   # need to add cell metadata (adata.obs / seuratObj@meta.data) separately
#   obs <- h5read(paste0(adata_path, '.h5seurat'), "/meta.data")
#   # format metadata as dataframe
#   meta <- data.frame(lapply(names(obs), function(x) {
#     if (length(obs[[x]])==2)
#       obs[[x]][['categories']][ifelse(obs[[x]][['codes']] >= 0, obs[[x]][['codes']] + 1, NA)]
#     else
#       as.numeric(obs[[x]])
#   }
#   ), row.names=Cells(seuratObj))
#   colnames(meta) <- names(obs)
# 
#   # add to seurat object
#   seuratObj <- AddMetaData(seuratObj,meta)
# 
#   saveRDS(seuratObj, paste0('Data/Derived/network/SCING_complete/supercells_seurat/', ct[i], ".RDS"))
# 
# }


