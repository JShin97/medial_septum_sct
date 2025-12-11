setwd("/u/project/xyang123/jshin/medial_septum_sct")

library(plyr)
library(dplyr)
library(tidyverse)
source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_Version_1.99.0.R")

celltype_networks <- "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/network/SCING_complete/saved_networks/final_edges/"
celltypes <- list.files(celltype_networks, full.names = FALSE, recursive = FALSE)
celltypes <- gsub(".csv.gz", "", celltypes)
effects <- c("Normative", "Gonad", "SexChrom", "Gonad:SexChrom", "addX", "addY1", "addY2")

run_KDA <- function(celltype, effect) {
  
  job.kda <- list()
  job.kda$label <- effect
  job.kda$folder <- paste0("Results/complete_analysis/KDA/DEG_SD/sexeffect_separate", "/", celltype, "/", effect) # parent folder for results
  # Input a network
  # columns: TAIL HEAD WEIGHT
  job.kda$netfile <- paste0("Results/complete_analysis/KDA/DEG_SD/input_files/network/", celltype, "_network.txt")
  # Tab delimited text file of gene set containing two columns: MODULE, NODE
  # Outputs from Module Merge script can be directly used
  job.kda$modfile <- paste0("Results/complete_analysis/KDA/DEG_SD/input_files/modfile/sexeffect_modules/", effect, "_", celltype, "_modfile.txt")
  
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
  
  job.kda <- kda2cytoscape(job.kda)
}

for(i in 7) {
  
  celltype <- celltypes[i]
  
  if(!dir.exists(paste0("Results/complete_analysis/KDA/DEG_SD/sexeffect_separate/", celltype))) {
    dir.create(paste0("Results/complete_analysis/KDA/DEG_SD/sexeffect_separate/", celltype))
  }
  
  try(run_KDA(celltype, "Normative"), silent = TRUE)
  try(run_KDA(celltype, "Gonad"), silent = TRUE)
  try(run_KDA(celltype, "SexChrom"), silent = TRUE)
  try(run_KDA(celltype, "Gonad:SexChrom"), silent = TRUE)
  try(run_KDA(celltype, "addX"), silent = TRUE)
  try(run_KDA(celltype, "addY1"), silent = TRUE)
  try(run_KDA(celltype, "addY2"), silent = TRUE)
} ##fix this so that it keeps going even after no KDs found


