##Script to run MSEA with mergeomics with Monty's old datasets + new datasets 
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
input_files <- c("modfile_autosomal.txt", "modfile_sexlinked.txt")
##############################################################################

# ## processing Monty's old datasets
# gwas_dir <- "/u/project/xyang123/jshin/resources/GWAS_genesets/"
# gwas_sets <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)
# mapping <- c("50.50") ##dist50 only
# 
# for(i in 67:length(gwas_sets)) {
# 
#   trait <- gwas_sets[i]
# 
#   for(j in 1:length(input_files)) {
# 
#     current_input <- input_files[j]
#     current_input_label <- current_input %>% gsub("modfile_", "", .) %>% gsub(".txt", "", .)
# 
#     try({job.ssea <- list()
#     job.ssea$label <- paste(gsub(".dist50", "", trait), mapping, sep = "_")
#     job.ssea$folder <- paste0("Results/complete_analysis/MSEA/DEG_SD/dist50/", current_input_label)
#     job.ssea$genfile <- paste0(gwas_dir, "/", trait, "/", mapping, ".g.txt")
#     job.ssea$locfile <- paste0(gwas_dir, "/", trait, "/", mapping, ".l.txt")
#     job.ssea$modfile <- paste0("Results/complete_analysis/MSEA/DEG_SD/input_files/", current_input)
#     #job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
#     job.ssea$permtype <- "gene" # for EWAS, set this to "marker"
#     job.ssea$nperm <- 5000
#     job.ssea$maxgenes <- 1000
#     job.ssea$maxoverlap <- 0.33 # for EWAS, set this to 1
#     job.ssea$trim <- 0.002 # default is 0.002, set to 0 for no trimming, users can try increasing to 0.005 to dampen inflation further
#     job.ssea <- ssea.start(job.ssea)
#     job.ssea <- ssea.prepare(job.ssea)
#     job.ssea <- ssea.control(job.ssea)
#     job.ssea <- ssea.analyze(job.ssea)
#     job.ssea <- ssea.finish(job.ssea)
#     }, silent = TRUE)
#   }
# }
# #############################################################################
# processing new GWAS datasets
gwas_dir <- "/u/project/xyang123/shared/datasets/neuro_gwas_2024/post_mdf/"
gwas_sets <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)

mapping <- list.files(paste0(gwas_dir, gwas_sets[1]), full.names = FALSE, recursive = FALSE)
mapping <- gsub(".g.txt|.m.txt", "", mapping) %>% unique()

for(i in 1:length(gwas_sets)) {

  trait <- gwas_sets[i]

  for(j in 1:length(input_files)) {

    current_input <- input_files[j]
    current_input_label <- current_input %>% gsub("modfile_", "", .) %>% gsub(".txt", "", .)

    try({job.ssea <- list()
    job.ssea$label <- paste(gsub(".dist50", "", trait), mapping, sep = "_")
    job.ssea$folder <- paste0("Results/complete_analysis/MSEA/DEG_SD/dist50/", current_input_label)
    job.ssea$genfile <- paste0(gwas_dir, trait, "/", mapping, ".g.txt")
    job.ssea$locfile <- paste0(gwas_dir, trait, "/", mapping, ".m.txt")
    job.ssea$modfile <- paste0("Results/complete_analysis/MSEA/DEG_SD/input_files/", current_input)
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
##############################################################################
##Processing updated gwas sets 
# completed_gwas_sets <- list.files("Results/celltype2_analysis/MSEA/dist50/logfc_abs0.1/msea/", pattern = ".pvalues.txt")
# completed_gwas_sets <- sapply(strsplit(completed_gwas_sets, "_"), function(x) paste(x[1:2], collapse = "_"))
# 
# gwas_dir1 <- "/u/scratch/j/jshin/GWAS_genesets/"
# gwas_sets1 <- list.dirs(gwas_dir1, full.names = FALSE, recursive = FALSE)
# gwas_sets1 <- gsub(".dist50", "", gwas_sets1)
# gwas_dir2 <- "/u/project/xyang123/shared/datasets/neuro_gwas_2024/post_mdf/"
# gwas_sets2 <- list.dirs(gwas_dir2, full.names = FALSE, recursive = FALSE)
# 
# gwas_sets_all <- c(gwas_sets1, gwas_sets2)
# gwas_sets <- gwas_sets_all[!gwas_sets_all %in% completed_gwas_sets]
# 
# mapping <- list.files(paste0(gwas_dir, gwas_sets[1]), full.names = FALSE, recursive = FALSE)
# mapping <- gsub(".g.txt|.m.txt", "", mapping) %>% unique()
# 
# for(i in 2:length(gwas_sets)) {
#   
#   trait <- gwas_sets[i]
#   
#   for(j in 1:length(input_files)) {
#     
#     current_input <- input_files[j]
#     current_input_label <- current_input %>% gsub("modfile_", "", .) %>% gsub(".txt", "", .)
#     
#     job.ssea <- list()
#     job.ssea$label <- paste(gsub(".dist50", "", trait), mapping, sep = "_")
#     job.ssea$folder <- paste0("Results/celltype2_analysis/MSEA/dist50/", current_input_label)
#     job.ssea$genfile <- paste0(gwas_dir, trait, "/", mapping, ".g.txt")
#     job.ssea$locfile <- paste0(gwas_dir, trait, "/", mapping, ".m.txt")
#     job.ssea$modfile <- paste0("Results/celltype2_analysis/MSEA/input_files/", current_input)
#     #job.ssea$inffile <- "../resources/genesets/kbr.info.txt"
#     job.ssea$permtype <- "gene" # for EWAS, set this to "marker"
#     job.ssea$nperm <- 5000
#     job.ssea$maxgenes <- 1000
#     job.ssea$maxoverlap <- 0.33 # for EWAS, set this to 1
#     job.ssea$trim <- 0.002 # default is 0.002, set to 0 for no trimming, users can try increasing to 0.005 to dampen inflation further
#     job.ssea <- ssea.start(job.ssea)
#     job.ssea <- ssea.prepare(job.ssea)
#     job.ssea <- ssea.control(job.ssea)
#     job.ssea <- ssea.analyze(job.ssea)
#     job.ssea <- ssea.finish(job.ssea)
#   }
# }


