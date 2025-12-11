##Genotype validation 

library(tidyverse)
library(limma)
library(edgeR)
library(data.table)
library(plyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(magrittr)
library(ggpubr)
library(ggplotify)
library(scran)
library(scater)
library(AnnotationHub)
library(DESeq2)
library(openxlsx)

seurat <- readRDS(file = "/u/scratch/j/jshin/cerebellum_sct/Saves/PostFilter/cellbender_clustered.RDS")

sample_gt_mapping <- openxlsx::read.xlsx("~/scratch/medial_septum_sct/Data/sample_genotype_mapping.xlsx") 
##fix Genotype metadata 
seurat_metadata <- seurat@meta.data %>% 
  as.data.frame() %>% 
  dplyr::select(SampleID, Genotype, Gonad, SexChrom, Trisomy) %>% 
  rownames_to_column(var = "cellID") %>% 
  left_join(sample_gt_mapping, by = c("SampleID", "Genotype")) %>% 
  column_to_rownames(var  = "cellID")
seurat@meta.data <- seurat_metadata
seurat$Sample2 <- gsub("M", "T", seurat$Sample2)
seurat$Sample2 <- gsub("F", "O", seurat$Sample2)
seurat$Sample2_original <- gsub("M", "T", seurat$Sample2_original)
seurat$Sample2_original <- gsub("F", "O", seurat$Sample2_original)

sce <- as.SingleCellExperiment(seurat)
sce.aggregate <- aggregateAcrossCells(sce, id = DataFrame(sce$Sample2))

#Genotype validation 
ah <- AnnotationHub()
ensdb <- ah[["AH73905"]]

chromosome <- mapIds(ensdb,
                     keys = rownames(sce.aggregate),
                     keytype = "SYMBOL",
                     column = "SEQNAME")
rowData(sce.aggregate)$chromosome <- chromosome
sce.aggregate <- sce.aggregate[!is.na(rowData(sce.aggregate)$chromosome), ]

#aggregate counts 
xist_counts <- as.numeric(counts(sce.aggregate["Xist", ]))
uty_counts <- as.numeric(counts(sce.aggregate["Uty", ]))
total_counts <- colSums(counts(sce.aggregate))

sexing_results <- data.frame(
  Sample = sce.aggregate$Sample2,
  Sample_original = sce.aggregate$Sample2_original, 
  Genotype = sce.aggregate$Genotype,
  Genotype_original = sce.aggregate$Genotype_original,
  TotalCounts = total_counts,
  XistCounts = xist_counts,
  UtyCounts = uty_counts
)

sexing_df <- sexing_results %>%
  pivot_longer(-c(Sample, Sample_original, Genotype, Genotype_original, TotalCounts)) %>%
  mutate(fraction = value / TotalCounts) %>%
  dplyr::select(-c(TotalCounts, value)) %>%
  pivot_wider(values_from = fraction, names_from = name)
sample_order <- sexing_df$Sample[order(sexing_df$Genotype)]
sexing_df <- sexing_df %>% dplyr::mutate(Sample = factor(Sample, levels = sample_order))
sample_order <- sexing_df$Sample_original[order(sexing_df$Genotype_original)]
sexing_df <- sexing_df %>% dplyr::mutate(Sample_original = factor(Sample_original, levels = sample_order))
sexing_df$Genotype <- gsub("F", "O", sexing_df$Genotype)
sexing_df$Genotype <- gsub("M", "T", sexing_df$Genotype)
sexing_df$Genotype_original <- gsub("F", "O", sexing_df$Genotype_original)
sexing_df$Genotype_original <- gsub("M", "T", sexing_df$Genotype_original)

plot1 <- sexing_df %>%
  ggplot(aes(x = XistCounts, y = Sample_original, fill = Genotype_original)) + geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  xlab("Fraction Xist / Total") +
  ylab("Sample") +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 10), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 15), 
        axis.title.y = element_text(size = 17), 
        axis.title.x = element_text(size = 17), 
        plot.title = element_text(size = 20), 
        legend.position = "none") + 
  ggtitle("Fraction Xist Count by Sample") 

plot2 <- sexing_df %>%
  ggplot(aes(x = UtyCounts, y = Sample_original, fill = Genotype_original)) + geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  xlab("Fraction Uty / Total") +
  ylab("Sample") + 
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 10), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 15), 
        axis.title.y = element_text(size = 17), 
        axis.title.x = element_text(size = 17), 
        plot.title = element_text(size = 20), 
        legend.position = "right") +
  ggtitle("Fraction Uty Count by Sample")

plot3 <- sexing_df %>%
  ggplot(aes(x = XistCounts, y = Sample, fill = Genotype)) + geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  xlab("Fraction Xist / Total") +
  ylab("Sample") +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 10), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 15), 
        axis.title.y = element_text(size = 17), 
        axis.title.x = element_text(size = 17), 
        plot.title = element_text(size = 20), 
        legend.position = "none") +
  ggtitle("Fraction Xist Count by Sample")

plot4 <- sexing_df %>%
  ggplot(aes(x = UtyCounts, y = Sample, fill = Genotype)) + geom_bar(position = "dodge", stat = "identity") +
  theme_bw() +
  xlab("Fraction Uty / Total") +
  ylab("Sample") +
  theme(axis.text.y = element_text(size = 15), 
        axis.text.x = element_text(size = 10), 
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 15), 
        axis.title.y = element_text(size = 17), 
        axis.title.x = element_text(size = 17), 
        plot.title = element_text(size = 20), 
        legend.position = "right") +
  ggtitle("Fraction Uty Count by Sample")
