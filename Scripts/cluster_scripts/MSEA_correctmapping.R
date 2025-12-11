##Script to run MSEA with mergeomics with re-processed gwas datasets 
##use older version of MSEA script 

setwd("/u/project/xyang123/jshin/medial_septum_sct")
library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_ZS_2020.r")

##############################################################################
## input files
input_dir <- "Results/complete_analysis/MSEA/DEG_SD_metacell/input_files"
#input_files <- list.files(input_dir, full.names = FALSE)
input_files <- c("modfile_logfc_abs0.1_top1000.txt")
#input_files <- c("modfile_neuron_top1000.txt")
##############################################################################

## processing re-processed gwas datasets
gwas_dir <- "/u/project/xyang123/rainyliu/temp/gwas_mdf/post_mdf"
gwas_sets <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)
mapping <- c("50.50") ##dist50 only

gwas_sets <- gwas_sets[1:(length(gwas_sets)-2)]

for(i in 99:length(gwas_sets)) {

  trait <- gwas_sets[i]

  for(j in 1:length(input_files)) {

    current_input <- input_files[j]
    current_input_label <- current_input %>% gsub("modfile_", "", .) %>% gsub(".txt", "", .)

    try({job.ssea <- list()
    job.ssea$label <- paste(gsub(".dist50", "", trait), mapping, sep = "_")
    job.ssea$folder <- paste0("Results/complete_analysis/MSEA/DEG_SD_metacell/dist50/", current_input_label)
    job.ssea$genfile <- paste0(gwas_dir, "/", trait, "/", mapping, ".g.txt")
    job.ssea$locfile <- paste0(gwas_dir, "/", trait, "/", mapping, ".l.txt")
    job.ssea$modfile <- paste0("Results/complete_analysis/MSEA/DEG_SD_metacell/input_files/", current_input)
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
    }, silent = TRUE)
  }
}
