##Script to run MSEA with mergeomics with Monty's old datasets 
##use older version of MSEA script 

setwd("/u/project/xyang123/jshin/medial_septum_sct")
library(plyr)
library(dplyr)
library(tidyverse)
source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_ZS_2020.r")

#If processing Monty's old datasets
gwas_dir <- "/u/home/j/jshin/scratch/GWAS_genesets/"
gwas_sets <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)
mapping <- list.files(paste0(gwas_dir, gwas_sets[1]), full.names = FALSE, recursive = FALSE)
mapping <- gsub(".g.txt|.l.txt", "", mapping) %>% unique()

for(i in 1:length(gwas_sets)) {
  
  trait <- gwas_sets[i]
  
  for(j in 1:length(mapping)) {
    
    current_mapping <- mapping[j]
    
    job.ssea <- list()
    job.ssea$label <- paste(gsub(".dist50", "", trait), current_mapping, sep = "_")
    job.ssea$folder <- "Results/MSEA/dist50_new/"
    job.ssea$genfile <- paste0(gwas_dir, "/", trait, "/", current_mapping, ".g.txt")
    job.ssea$locfile <- paste0(gwas_dir, "/", trait, "/", current_mapping, ".l.txt")
    job.ssea$modfile <- "Results/MSEA/input_files/modfile_logfc0.1.txt"
    #job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
    job.ssea$permtype <- "gene" # for EWAS, set this to "marker"
    job.ssea$nperm <- 5000
    job.ssea$maxgenes <- 1000
    job.ssea$maxoverlap <- 0.33 # for EWAS, set this to 1
    job.ssea$trim <- 0.002 # default is 0.002, set to 0 for no trimming, users can try increasing to 0.005 to dampen inflation further
    job.ssea <- ssea.start(job.ssea)
    job.ssea <- ssea.prepare(job.ssea)
    job.ssea <- ssea.control(job.ssea)
    job.ssea <- ssea.analyze(job.ssea)
    job.ssea <- ssea.finish(job.ssea)
    
  }
}