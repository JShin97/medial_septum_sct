##Script to run kDA analysis on cell-type specific networks

setwd("/u/project/xyang123/jshin/medial_septum_sct")

library(plyr)
library(dplyr)
library(tidyverse)
source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_Version_1.99.0.R")

# celltype_networks <- "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/network/SCING_complete/saved_networks/final_edges/"
# celltypes <- list.files(celltype_networks, full.names = FALSE, recursive = FALSE)
# celltypes <- gsub(".csv.gz", "", celltypes)
celltypes <- c("Astrocyte", "Ependymal", "GABA-Chol", "GABA", "Glut", "Inh_IMN", "Microglia", "Oligo", "OPC", "Vascular")
effects <- c("Normative", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2")

##Step 2 : run wKDA analysis -- all effects together
# modules <- list.files("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/modfile/celltype_modules", full.names = FALSE, recursive = FALSE)
# modules <- gsub("_modfile.txt", "", modules)
# for(i in 5) {
#   
#   module <- modules[i]
#   # celltype <- str_split(module, pattern = "_", n = 2)[[1]][1]
#   # effect <- str_split(module, pattern = "_", n = 2)[[1]][2]
#   
#   job.kda <- list()
#   job.kda$label <- module
#   job.kda$folder <- paste0("Results/complete_analysis/KDA/DEG_SD_metacell/", module) # parent folder for results
#   # Input a network
#   # columns: TAIL HEAD WEIGHT
#   job.kda$netfile <- paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/network/", module, "_network.txt")
#   # Tab delimited text file of gene set containing two columns: MODULE, NODE
#   # Outputs from Module Merge script can be directly used
#   job.kda$modfile <- paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/modfile/celltype_modules/", module, "_modfile.txt")
#   
#   # 0.0-1.0 - 0.0 means not factoring in edge weights, 0.5 means partial influence,
#   # and 1.0 means full influence
#   job.kda$edgefactor<-0 #or 0 (less restrictions than 1)
#   # The searching depth for the KDA
#   job.kda$depth<-1 #can try 2 
#   # 0 means we do not consider the directions of the regulatory interactions
#   # while 1 is opposite.
#   job.kda$direction <- 0 
#   job.kda$nperm<-10000 ##10000 is better but takes longer 
#   
#   ## Let's run KDA!
#   job.kda <- kda.configure(job.kda)
#   job.kda <- kda.start(job.kda)
#   job.kda <- kda.prepare(job.kda)
#   job.kda <- kda.analyze(job.kda)
#   job.kda <- kda.finish(job.kda)
#   
#   job.kda <- kda2cytoscape(job.kda)
#   
# }

##Step 2 : run wkDA analysis -- each effect separately 
for(i in 1:length(celltypes)) {
  celltype <- celltypes[i]
  for(j in 1:length(effects)) {
    
    try({effect <- effects[j]
    module <- paste0(effect, "_", celltype)
    
    job.kda <- list()
    job.kda$label <- celltype
    job.kda$folder <- paste0("Results/complete_analysis/KDA/DEG_SD_metacell/sexeffect_results/", celltype, "/", effect) # parent folder for results -- make celltype folder manually
    # Input a network
    # columns: TAIL HEAD WEIGHT
    job.kda$netfile <- paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/network/", celltype, "_network.txt")
    # Tab delimited text file of gene set containing two columns: MODULE, NODE
    # Outputs from Module Merge script can be directly used
    job.kda$modfile <- paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/modfile/sexeffect_modules/", module, "_modfile.txt")
    
    # 0.0-1.0 - 0.0 means not factoring in edge weights, 0.5 means partial influence,
    # and 1.0 means full influence
    job.kda$edgefactor<-0 #or 0 (less restrictions than 1)
    # The searching depth for the KDA
    job.kda$depth<-1 #can try 2
    # 0 means we do not consider the directions of the regulatory interactions
    # while 1 is opposite.
    job.kda$direction <- 0
    job.kda$nperm<-10000 ##10000 is better but takes longer
    
    ## Let's run KDA!
    job.kda <- kda.configure(job.kda)
    job.kda <- kda.start(job.kda)
    job.kda <- kda.prepare(job.kda)
    job.kda <- kda.analyze(job.kda)
    job.kda <- kda.finish(job.kda)
    
    job.kda <- kda2cytoscape(job.kda)}, silent = TRUE)
  }
}
