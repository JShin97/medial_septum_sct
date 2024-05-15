setwd("/u/project/xyang123/jshin/medial_sepum_sct")

library(tidyverse)
library(dplyr)
library(stringr)

original <- readRDS("Saves/Caden/For_Figures.rds")
sample_metadata <- original@meta.data %>% select(orig.ident, Genotype, Gonad, SexChrom, Trisomy) %>% unique()
rownames(sample_metadata) <- NULL

sample_mapping <- read_csv("Data/Raw/sample_mapping.csv") %>% select(Sample_ID, Sample_Name)
sample_mapping$Sample_ID <- gsub("XY-", "", sample_mapping$Sample_ID)
sample_mapping$Sample_Name <- gsub("XY-", "", sample_mapping$Sample_Name)
sample_mapping$Sample_ID <- str_sub(sample_mapping$Sample_ID,1, -3)
sample_mapping$Sample_Name <- str_sub(sample_mapping$Sample_Name,1, -3)
sample_mapping <- sample_mapping %>% unique()

technical_metadata <- readxl::read_xlsx("Data/Raw/technical_metadata.xlsx", sheet = "10x-processing") %>% select(`10X run`, `Lib prep`, `Sample ID...5`, `Genotype...6`)
colnames(technical_metadata) <- c("Sequencing_Date", "Library_Prep", "Sample_Name", "Genotype")
technical_metadata$Sample_Name <- paste0("TriSeptum-", technical_metadata$Sample_Name)
technical_metadata <- left_join(technical_metadata, sample_mapping, by = c("Sample_Name")) %>% relocate(Sample_Name)
technical_metadata$Sequencing_Date <- factor(technical_metadata$Sequencing_Date)
technical_metadata$Library_Prep <- factor(technical_metadata$Library_Prep)

save(sample_metadata, sample_mapping, technical_metadata, file = "Saves/metadata.Rda") 
