##Prep for MapMyCells 

library(Seurat)
library(data.table)
library(rhdf5)
library(anndata)
library(magrittr)

#########################################################
##Part 1 : Prep input (run in jupyter notebook instead)
#########################################################
##Rows are expected to be cells and columns are expected to be genes
seuratObject <- readRDS("/u/scratch/j/jshin/medial_septum_sct/Saves/For_Figures.rds")
counts <- seuratObject@assays$RNA@counts
counts <- counts %>% t()

genes   <- colnames(counts)
samples <- rownames(counts)
sparse_counts <- as(counts, "dgCMatrix")  # Convert to sparse matrix, if not already
write.csv(sparse_counts, file = "/u/scratch/j/jshin/medial_septum_sct/Data/Derived/counts.csv", row.names = TRUE, col.names = TRUE)

countAD <- anndata::AnnData(X   = sparse_counts,   # Create the anndata object
                   var = data.frame(genes=genes,row.names=genes),
                   obs = data.frame(samples=samples,row.names=samples))
write_h5ad(countAD, "counts.h5ad") # Write it out as h5ad

#########################################################
##Part 2 : Run MapMyCells (website)
#########################################################

#########################################################
##Part 3 : Process mapping results 
#########################################################
##Load Seurat object 
seuratObject <- readRDS("/u/scratch/j/jshin/medial_septum_sct/Saves/For_Figures_clustered.rds")
metadata <- seuratObject@meta.data ##rownames = cells
cellids <- rownames(metadata)

##Load mapping results
mapping <- read.csv("Data/Derived/MapMyCells/medial_septum_for_MapMyCells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1724376165222.csv",
                    comment.char="#")
head(data.frame(mapping))
data.frame(Cell_counts=head(sort(table(mapping$class_name),decreasing=T),8))
sum(t(t(head(sort(table(mapping$class_name),decreasing=T),8))))/length(mapping$class_name)
data.frame(Cell_counts = sort(table(mapping$class_name), decreasing = T)) 
data.frame(Cell_counts = sort(table(mapping$subclass_name), decreasing = T)) 

cellcounts <- table(mapping$class_name) %>% as.data.frame()
subcellcounts <- table(mapping$subclass_name) %>% as.data.frame()

#Assign cells with low probability as "unknown"
#mapping[mapping$cluster_bootstrapping_probability < 0.90, "class_new"] <- "Unknown"
#mapping[mapping$cluster_bootstrapping_probability < 0.90, "subclass_new"] <- "Unknown"

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
mapping[grepl("Neuronal", mapping$NT), "NT"] <- "Misc_NT"
#mapping[grepl("Other", mapping$class_new), "NT"] <- "Other"
#mapping[grepl("Other", mapping$subclass_new), "NT"] <- "Other"

#Simplify by broad class 
mapping$class_simplify <- mapping$NT
mapping[mapping$NT == "Non-neuronal", "class_simplify"] <- mapping[mapping$NT == "Non-neuronal", "class_new"]
mapping$class_simplify <- plyr::mapvalues(mapping$class_simplify,
                                          from = c("30 Astro-Epen", "31 OPC-Oligo", "33 Vascular", "34 Immune"),
                                          to = c("Astro-Epen", "OPC-Oligo", "Vascular", "Immune")
                                          )
mapping[mapping$class_new == "Other", "class_simplify"] <- "Other"

#Simplify by subclass
mapping$subclass_simplify <- mapping$NT
mapping[mapping$NT == "Non-neuronal", "subclass_simplify"] <- mapping[mapping$NT == "Non-neuronal", "subclass_new"]
mapping$subclass_simplify <- plyr::mapvalues(mapping$subclass_simplify,
                                             from = c("318 Astro-NT NN", "319 Astro-TE NN", "321 Astroependymal NN", "323 Ependymal NN", "325 CHOR NN",
                                                      "326 OPC NN", "327 Oligo NN", "330 VLMC NN", "331 Peri NN", "333 Endo NN", "334 Microglia NN"),
                                             to = c("Astro-NT", "Astro-TE", "Astroependymal", "Ependymal", "Choroid", "OPC", "Oligo", "VLMC", "Pericyte", 
                                                    "Endothelial", "Microglia"))
mapping[mapping$subclass_new == "Other", "subclass_simplify"] <- "Other"
#mapping[grepl("NN", mapping$subclass_simplify), "subclass_simplify"] %>% stringr::str_split_fixed(., " ", n = 3)

#Assign cells with low probability as "unknown"
#mapping[mapping$cluster_bootstrapping_probability < 0.50, "class_simplify"] <- "Unknown"
#mapping[mapping$cluster_bootstrapping_probability < 0.50, "subclass_simplify"] <- "Unknown"

#Put row.names as data colnames and the order to match the data
rownames(mapping) <- mapping$cell_id
mapping <- mapping[cellids,]
#metadata <- cbind(metadata, mapping)

#Add new metadata columns from mapping 
seuratObject <- AddMetaData(seuratObject, 
                            metadata = mapping,
                            col.name = colnames(mapping))
#seuratObject_subset <- subset(seuratObject, class_new %in% c("other", "unknown"), invert = TRUE)

# celltype_results <- as.data.frame.matrix(table(seuratObject$cell.type, seuratObject$class_new))
# write.csv(celltype_results, file = "Results/MapMyCells/celltype_old_new.csv", row.names = TRUE)
# subcelltype_results <- as.data.frame.matrix(table(seuratObject$subtype, seuratObject$class_new))
# write.csv(subcelltype_results, file = "Results/MapMyCells/subcelltype_old_new.csv", row.names = TRUE)

#Visualization 
DimPlot(seuratObject, reduction = "umap", group.by = "neuronal")
DimPlot(seuratObject, reduction = "umap", group.by = "NT")
DimPlot(seuratObject, reduction = "umap", group.by = "class_simplify")
DimPlot(seuratObject, reduction = "umap", group.by = "subclass_simplify")

Idents(seuratObject) <- "class_simplify"
other_cells <- WhichCells(seuratObject, idents = c("Unknown"))
DimPlot(seuratObject, reduction = "umap", cells.highlight = other_cells, raster = FALSE)

DimPlot(seuratObject_subset, reduction = "umap", group.by = "class_new")
Idents(seuratObject_subset) <- "class_new"
other_cells <- WhichCells(seuratObject_subset, idents = c("34 Immune"))
DimPlot(seuratObject_subset, reduction = "umap", cells.highlight = other_cells, raster = FALSE)

p1 <- DimPlot(seuratObject, reduction = "umap", group.by = "class_new")
p2 <- DimPlot(seuratObject, reduction = "umap", group.by = "subclass_new")
p3 <- DimPlot(seuratObject, reduction = "umap", group.by = "cell.type")
p4 <- DimPlot(seuratObject, reduction = "umap", group.by = "subtype")

p3 + p4
p1 + p2

png("Plots/MapMyCells/new_class.png", res = 600, width = 5, height = 5)
print(p1)
dev.off()

##Save
saveRDS(seuratObject, file = "/u/scratch/j/jshin/medial_septum_sct/Saves/For_Figures_clustered2.rds")

#########################################################
##Part 4 : Doublet detection (DF.Rmd)
#########################################################

#########################################################
##Part 5 : Create final seurat object
#########################################################
seuratObject <- readRDS(file = "/u/scratch/j/jshin/medial_septum_sct/Saves/For_Figures_clustered2.rds")
seuratObject <- subset(seuratObject, doublet_cluster == "FALSE")

Idents(seuratObject) <- "subclass_simplify"
cells <- WhichCells(seuratObject, idents = c("Ependymal"))
DimPlot(seuratObject, reduction = "umap", cells.highlight = cells, raster = FALSE)

metadata$cell.type2 <- metadata$subtype
metadata[metadata$subtype == "Unknown", "cell.type2"] <- "Astrocyte"
metadata[metadata$subtype %in% c("Vglut1", "Vglut2"), "cell.type2"] <- "Glut"
metadata[metadata$subtype %in% c("GABAergic"), "cell.type2"] <- "GABA"
metadata[metadata$subtype %in% c("Cholinergic") & metadata$subclass_simplify %in% c("GABA-Chol"), "cell.type2"] <- "GABA-Chol"
metadata[metadata$subtype %in% c("Cholinergic") & !metadata$subclass_simplify %in% c("GABA-Chol"), "cell.type2"] <- "GABA"
metadata[metadata$subtype == "Intermediate Progenitor Cell", "cell.type2"] <- "Inh_IMN"
metadata[metadata$subtype == "Oligodendrocyte Progenitor Cell", "cell.type2"] <- "OPC"
metadata[metadata$subtype == "Myelinating Oligodendrocyte", "cell.type2"] <- "Oligo"
metadata[metadata$subtype == "Endothelial", "cell.type2"] <- "Vascular"
metadata[metadata$subtype == "Endothelial" & metadata$subclass_simplify %in% c("GABA"), "cell.type2"] <- "GABA"
metadata[metadata$subtype == "Microglia" & metadata$subclass_simplify %in% c("GABA"), "cell.type2"] <- "GABA"
metadata[metadata$subtype == "Oligodendrocyte Progenitor Cell" & metadata$subclass_simplify == c("GABA"), "cell.type2"] <- "GABA"
metadata[metadata$subtype == "Oligodendrocyte Progenitor Cell" & metadata$subclass_simplify == c("Oligo"), "cell.type2"] <- "Oligo"
metadata[metadata$subtype == "Myelinating Oligodendrocyte" & metadata$subclass_simplify == c("GABA"), "cell.type2"] <- "GABA"

seuratObject$cell.type2 <- metadata$cell.type2
DimPlot(seuratObject, group.by = "cell.type2")

saveRDS(seuratObject, file = "/u/scratch/j/jshin/medial_septum_sct/Saves/For_Figures_new.rds")







