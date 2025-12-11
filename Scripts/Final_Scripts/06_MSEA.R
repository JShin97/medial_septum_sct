##Script to create input files for MSEA and visualize results 

library(plyr)
library(dplyr)
library(tidyverse)
library(Seurat)
library(biomaRt)
library(ComplexHeatmap)
library(ggpubr)
library(RColorBrewer)
library(stringr)

source("./source/color_palette.R")
source("./source/Mergeomics_ZS_2020.r")

########################################################################
## Part 1 : MSEA input prep
########################################################################

# Function for conversion of mouse to human genes 
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep="\t")

convert_mouse_to_human <- function(gene_list){
  
  output = list()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% dplyr::filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% dplyr::filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      output[[gene]] <- data.frame(rep(gene, length(human_genes)), human_genes)
      colnames(output[[gene]]) <- c("mouse_gene", "human_gene")
    }
  }
  return (output)
}

##prepare modfile 
load(file = "Results/complete_analysis/DEG_processed/DEG_SD_metacell/DEG_df_RNA_SD.rda")

DEG_df <- DEG_df %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.1)

gene_annotation <- readRDS("Gene_annotation/gene_annotation.Rds")
gene_annotation2 <- gene_annotation[!duplicated(gene_annotation$mgi_symbol), c(1:3, 7:9)]
gene_universe <- read.table("Gene_annotation/gene_universe.txt") %>% unlist(use.names = FALSE)

DEG_df <- left_join(DEG_df, gene_annotation2, by = c("gene"="mgi_symbol"))
DEG_df$chromosome_name <- factor(DEG_df$chromosome_name, levels = c((1:19),"X","Y","MT"))

mouse_genes <- unique(DEG_df$gene)
orthologs <- convert_mouse_to_human(mouse_genes)
ortholog_df <- do.call("rbind", orthologs)

MSEA_input <- DEG_df %>% 
  mutate(module = paste(effect, cell_type, sep = "_")) %>%
  left_join(ortholog_df, by = c("gene" = "mouse_gene")) %>%
  group_by(module) %>%
  dplyr::arrange(adj.P.Val, .by_group = TRUE) %>%
  slice_head(n = 1000) %>%
  dplyr::select(module, human_gene)
colnames(MSEA_input) <- c("MODULE", "GENE")
MSEA_input <- MSEA_input[!is.na(MSEA_input$GENE), ]

write.table(MSEA_input, file = "Results/complete_analysis/MSEA/DEG_SD_metacell/input_files/modfile_logfc_abs0.1_top1000.txt", sep = "\t", quote = FALSE)

########################################################################
## Part 2 : run MSEA 
########################################################################

## input files
input_dir <- "Results/complete_analysis/MSEA/DEG_SD_metacell/input_files"
input_files <- c("modfile_logfc_abs0.1_top1000.txt")

gwas_dir <- "/u/project/xyang123/rainyliu/temp/gwas_mdf/post_mdf"
gwas_sets <- list.dirs(gwas_dir, full.names = FALSE, recursive = FALSE)
mapping <- c("50.50") 

gwas_sets <- gwas_sets[1:(length(gwas_sets)-2)]

for(i in 1:length(gwas_sets)) {
  
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

########################################################################
## Part 3 : Visualize results
########################################################################

##get traits  
trait_dir <- "Results/complete_analysis/MSEA/DEG_SD_metacell/dist50/logfc_abs0.1_top1000/msea/"
alltraits <- list.files(trait_dir, full.names = FALSE, recursive = FALSE)
alltraits <- sapply(alltraits, function(x) sub("_[^_]*$", "", x)) %>% unique()

cutoff <- "logfc_abs0.1_top1000"

results <- list()
for(i in 1:length(alltraits)) {
  
  trait <- alltraits[i]
  
  trait_result <- read.delim(Sys.glob(paste0(data_dir, "/", cutoff, "/msea/", trait, "_",  "*.pvalues.txt")))
  trait_result <- trait_result[, c(1,2)] ##2 for raw p value, 3 for FDR 
  colnames(trait_result) <- c("MODULE", as.character(trait))
  results[[trait]] <- trait_result
}

result_merged <- Reduce(function(x,y) merge(x, y, by = "MODULE"), results)
# rename_map <- c(
#   "Gonad:SexChrom" = "GExSCE",
#   "SexChrom"       = "SCE",
#   "Gonad"          = "GE",
#   "addX"           = "XCD",
#   "addY1"          = "YCD1",
#   "addY2"          = "YCD2", 
#   "Normative"      = "Typical")
# result_merged$MODULE <- str_replace_all(result_merged$MODULE, rename_map)

## Rename studies to more accurate naming 
result_merged <- result_merged %>% 
  dplyr::rename(CARDIOGRAMC4D = "Aragam_Coronary_artery_disease", 
                Davies_MemoryPerf = "Davies_Learning", 
                PAGE_FastingGlucoseBMIadj = "Downie_BMI_adjusted_fasted_blood_glucose", 
                Graham_HDL = "Graham_HDL_cholesterol_measurement", 
                Graham_LDL = "Graham_LDL_cholesterol_measurement", 
                ILAE_GenEpilepsy = "ILAE_Epilepsy", 
                AABCG_BreastCancer = "Jia_Breast_cancer", 
                Levin_HeartFailure = "Levin_Heart_failure", 
                Nagel_DepressedAffect = "Nagel_Depression", 
                PGC_ASD_SCZ = "PGC_ASD", 
                ReproGen_AgeAtMenopause = "ReproGen_AgeAtMenarche", 
                Sakaue_PsoriasisVulgaris = "Sakaue_Psoriasis_vulgaris", 
                Sakaue_SubstanceDependence = "Sakaue_Substance_Dependence", 
                Savage_Intelligence = "Savage_Cognition", 
                Shen_PTSD = "Shen_Anxiety", 
                Soo_EpisodicMemory = "Soo_Memory", 
                Wightman_AD_LateOnset = "Wightman_Late-onset_Alzheimer", 
                Woodbury_EmotionRecognition = "Woodbury_FacialRecognition", 
                Zhou_NicotineDependence = "Zhou_Nicotine_Dependence", 
                GOT2D_T2D = "DIAGRAM_T2D") %>% 
  dplyr::select(-UKBB_CAD)

#calculate adjusted p-values across diseases
p_value_columns <- result_merged[, 2:ncol(result_merged)]
pooled_p_values <- as.vector(as.matrix(p_value_columns))

# apply the Benjamini-Hochberg (BH) correction
# The total number of tests (N) is automatically inferred by the length of the vector.
q_values <- p.adjust(pooled_p_values, method = "BH")

# reshape the resulting q-values back into the original matrix structure
q_value_matrix <- matrix(q_values, 
                         nrow = nrow(p_value_columns), 
                         ncol = ncol(p_value_columns))

result_merged_unadjusted <- result_merged
result_merged[, 2:ncol(result_merged)] <- q_value_matrix

rownames(result_merged) <- result_merged$MODULE
result_merged$MODULE <- NULL
result_merged <- as.matrix(result_merged) %>% t()
result_merged <- result_merged[, !colnames(result_merged) %in% c("_ctrlA", "_ctrlB")]
result_merged <- result_merged[, !grepl("Misc", colnames(result_merged))]
colnames(result_merged) <- gsub("GABA-Chol", "Chol", colnames(result_merged))
colnames(result_merged) <- gsub("Inh_IMN", "Inh-IMN", colnames(result_merged))
result_merged <- -log10(pmax(result_merged, 1e-05))

# Define the order patterns and cell types
order_patterns <- c("Typical", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2")
cell_types <- c("GABA", "Chol", "Glut", "Oligo", "OPC", "Astrocyte", "Ependymal", "Microglia", "Vascular", "Inh-IMN")

vec <- colnames(result_merged)
split_vec <- strsplit(vec, "_")
# Extract the first and second variables
first_vars <- sapply(split_vec, function(x) x[1])
second_vars <- sapply(split_vec, function(x) paste(x[-1], collapse = "_"))

# Convert first variables to factors for ordering
first_vars <- factor(first_vars, levels = order_patterns)

# Convert second variables to factors for ordering
second_vars <- factor(second_vars, levels = cell_types)

# Create a dataframe for sorting
df <- data.frame(
  Original = vec,
  FirstVar = first_vars,
  SecondVar = second_vars,
  stringsAsFactors = FALSE
)

# Sort the dataframe based on FirstVar and SecondVar
df_sorted <- df[order(df$FirstVar, df$SecondVar), ]

# Extract the sorted vector
vec_sorted <- df_sorted$Original
df_reordered <- result_merged[, vec_sorted]
result_merged <- df_reordered

# Split column names to get cell type and effect
split_columns <- strsplit(colnames(df_reordered), "_")
effects <- sapply(split_columns, function(x) x[1])
cell_types <- sapply(split_columns, function(x) paste(x[-1], collapse = "_"))

# Create color maps for cell types and effects
cell_type_colors <- cell_hex_colors[!names(cell_hex_colors) %in% c("Misc-NT")]
effect_colors <- c("Typical" = "#FC8D62",
                   "GE" = "#E78AC3",
                   "SCE" = "#A6D854",
                   "GExSCE" = "#FFD92F",
                   "XCD" = "#78b4ff",
                   "YCD1" = "#66C2A5",
                   "YCD2" = "#b4b4ff" )

colAnnotation <- HeatmapAnnotation(
  CellType = cell_types,
  Effect = effects,
  col = list(
    CellType = cell_type_colors,
    Effect = effect_colors
  ),
  annotation_name_gp = grid::gpar(fontsize = 18),
  annotation_legend_param = list(
    title_gp = grid::gpar(fontsize = 18),      # Legend title font size
    labels_gp = grid::gpar(fontsize = 18))
  #show_legend = FALSE
)

plot <- Heatmap(result_merged, name = "-log10(pvalue)", col = c("white", "red"), 
                top_annotation = colAnnotation,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                show_column_names = TRUE,
                show_row_names = TRUE,
                row_names_gp = grid::gpar(fontsize = 20),
                column_names_gp = grid::gpar(fontsize = 20),
                heatmap_legend_param = list(
                  title_gp = grid::gpar(fontsize = 18),      # Legend title font size
                  labels_gp = grid::gpar(fontsize = 18)    # Legend labels font size
                ),
                row_names_max_width = unit(100, "cm"),
                column_names_max_height = unit(100, "cm"),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(result_merged[i, j] > -log10(0.001)) {
                    grid.text("***", x, y)
                  } else if(result_merged[i, j] > -log10(0.01)) {
                    grid.text("**", x, y)
                  } else if(result_merged[i, j] > -log10(0.05)) {
                    grid.text("*", x, y)
                  }
                })

plot <- draw(plot, merge_legends = TRUE, annotation_legend_side = "left", heatmap_legend_side = "left")
