##Script to run meta-MSEA 

setwd("/u/project/xyang123/jshin/medial_septum_sct")
library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)

#source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_ZS_2020.r")
source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_Version_1.99.0.R")
source("/u/project/xyang123/jshin/resources/mergeomics/Mergeomics_utils.R")

# prep inputs
gwas_dir <- "/u/project/xyang123/rainyliu/temp/gwas_mdf/post_mdf/"
gwas_sets_1 <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)
mapping <- c("50.50") ##dist50 only

marker_associations_1 = paste0(gwas_dir, gwas_sets_1, "/", mapping, ".l.txt")
marker_mappings_1 = paste0(gwas_dir, gwas_sets_1, "/", mapping, ".g.txt")

gwas_dir <- "/u/project/xyang123/shared/datasets/neuro_gwas_2024/post_mdf/"
gwas_sets_2 <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)
mapping <- list.files(paste0(gwas_dir, gwas_sets_2[1]), full.names = FALSE, recursive = FALSE)
mapping <- gsub(".g.txt|.m.txt", "", mapping) %>% unique()

marker_associations_2 = paste0(gwas_dir, gwas_sets_2, "/", mapping, ".m.txt")
marker_mappings_2 = paste0(gwas_dir, gwas_sets_2, "/", mapping, ".g.txt")

marker_associations = c(marker_associations_1, marker_associations_2)
marker_mappings = c(marker_mappings_1, marker_mappings_2)
gwas_sets = c(gwas_sets_1, gwas_sets_2) %>% gsub(".dist50", "", .)

# Option 2:
# Running Meta-MSEA on finished MSEA runs returned by runMSEA()
modfile <- "Results/complete_analysis/MSEA/DEG_SD/input_files/modfile_logfc_abs0.1_top1000.txt"

# msea_job1 <- runMSEA(#job = msea_job, 
#                      marker_set = modfile, 
#                      marker_associations = marker_associations[1], 
#                      marker_mapping = marker_mappings[1],
#                      output_dir="Results/complete_analysis/MSEA/DEG_SD/dist50/neuron_top1000_meta/", 
#                      label=gwas_sets[1])
# msea_job2 <- runMSEA(marker_set = modfile, 
#                      marker_associations = marker_associations[2], 
#                      marker_mapping = marker_mappings[2], 
#                      output_dir="Results/complete_analysis/MSEA/DEG_SD/dist50/neuron_top1000_meta/", 
#                      label=gwas_sets[2])
# 
# jobs <- list(msea_job1, msea_job2)

jobs <- lapply(seq_along(marker_associations), function(i) {

  job <- runMSEA(marker_set = modfile,
                 marker_associations = marker_associations[i],
                 marker_mapping = marker_mappings[i], 
                 output_dir = "Results/complete_analysis/MSEA/DEG_SD/dist50/neuron_top1000_meta/", 
                 label=gwas_sets[i])

  return(job)

})
job.meta <- ssea.meta(jobs, label = "metaMSEA", folder = "Results/complete_analysis/MSEA/DEG_SD/dist50/neuron_top1000_meta")
job.meta.kda <- list()
job.meta.kda$metamsea <- job.meta
job.meta.kda$modfile <- modfile

# Next, run KDA with completed Meta-MSEA job
# kda_job <- runKDA(job=job.meta.kda, 
#                   MSEA_fdr_cutoff = c(0.5,0.5,0.5),
#                   merge_modules=TRUE, 
#                   network="networks/bayesian.hs.blood.txt") 
