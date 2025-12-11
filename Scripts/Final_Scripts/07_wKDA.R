library(plyr)
library(dplyr)
library(tidyverse)
library(gtools)

source("./source/Mergeomics_Version_1.99.0.R")
source("./source/color_palette.R")

######################################################################
## Step 1 : Read in constructed SCING networks 
######################################################################
##Step 1 : Read in input variables
celltype_networks <- "/u/home/j/jshin/scratch/medial_septum_sct/Data/Derived/network/SCING_complete/saved_networks/final_edges/"
celltypes <- list.files(celltype_networks, full.names = FALSE, recursive = FALSE)
celltypes <- gsub(".csv.gz", "", celltypes)

######################################################################
## Step 2 : Prepare network file to be readable by mergeomics
######################################################################
##Output = celltype-specific network and module file 

for(i in 1:length(celltypes)) {
  
  celltype <- celltypes[i]
  
  network <- read.csv2(paste0(celltype_networks, celltype, ".csv.gz"), sep = ",", row.names = 1)
  colnames(network) <- c("TAIL", "HEAD", "WEIGHT")

  write.table(network, file = paste0("Results/complete_analysis/KDA/DEG_SD/input_files/network/", celltype, "_network.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

######################################################################
## Step 2.5 : Prepare module file to be readable by mergeomics 
######################################################################
##prepare modfile 
##Limma results 
load(file = "Results/complete_analysis/DEG_processed/DEG_SD_metacell/DEG_df_RNA_SD.rda")
DEG_df <- DEG_df %>% dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 0.1)
DEG_list <- split(DEG_df, f = DEG_df$effect)

wKDA_input <- DEG_df %>% 
  mutate(module = paste(effect, cell_type, sep = "_")) %>%
  dplyr::select(module, gene)
colnames(wKDA_input) <- c("MODULE", "NODE")
wKDA_input <- wKDA_input[!is.na(wKDA_input$NODE), ]

##filter out mitochondrial genes 
wKDA_input <- wKDA_input[!grepl(pattern = "mt-", wKDA_input$NODE), ]
celltypes <- unique(DEG_df$cell_type)

## Prepare celltype and sex-effect specific modfile 
effects <- unique(DEG_df$effect)
for(i in 1:length(celltypes)) {
  
  celltype <- celltypes[i]
  for(j in 1:length(effects)) {
    
    effect <- effects[j]
    module <- paste(effect, celltype, sep = "_")
    module_modfile <- wKDA_input %>% dplyr::filter(MODULE %in% module)
    
    write.table(module_modfile, file = paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/modfile/sexeffect_modules/", module, "_modfile.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
}

######################################################################
##Step 3 : run wkDA analysis considering each effect separately 
######################################################################
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

##################################################################################
##Step 5 : Create new network and node file based on KDA results to export to Cytoscape
##################################################################################
effects <- c("Typical", "GE", "SCE", "GExSCE", "XCD", "YCD1", "YCD2")
celltypes <- c("Astrocyte", "Ependymal", "GABA-Chol", "GABA", "Glut", "Inh_IMN", "Microglia", "Oligo", "OPC", "Vascular")

celltype <- "Glut"
celltype_modfile <- read.delim(file = paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/modfile/celltype_modules/", celltype, "_modfile.txt"))
celltype_modfile <- celltype_modfile %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = MODULE, values_from = value, values_fill = 0)

kd_list <- list()
effect_results <- list.dirs(paste0("Results/complete_analysis/KDA/DEG_SD_metacell/sexeffect_separate/", celltype), full.names = FALSE)
effect_results <- effect_results[grepl("cytoscape", effect_results)]
effect_results <- gsub("/cytoscape", "", effect_results)
for(j in 1:length(effect_results)) {
  
  effect <- effect_results[j]
  kds <- read.delim(paste0("Results/complete_analysis/KDA/DEG_SD_metacell/sexeffect_separate/", celltype, "/", effect, "/kda/", celltype, ".pvalues.txt"))
  colnames(kds) <- c("MODULE", "NODE", "Pval", "FDR", "FOLD")
  kds <- kds %>% 
    dplyr::filter(FDR < 0.05) %>% 
    arrange(Pval)
  top_kds <- kds %>% 
    mutate(effect =  gsub(paste0("_", celltype), "", MODULE)) %>%
    dplyr::select(effect, NODE, FDR) %>% 
    #mutate(value = 1) %>%
    pivot_wider(names_from = effect, values_from = FDR, values_fill = 0) ##fill with Pval and not value 
  slice_head(n = 5)
  kd_list[[effect]] <- top_kds
  
}
kd_list <- kd_list[effects]
names(kd_list) <- effects

# Replace NULL elements with an empty dataframe 
kd_list <- Map(function(df, name) {
  if (is.null(df)) {
    data.frame(NODE = character(), setNames(list(double()), name))
  } else {
    df
  }
}, kd_list, names(kd_list))

# Function to merge dataframes in list by 'genes'
kd_df <- Reduce(function(x, y) full_join(x, y, by = "NODE"), kd_list)
kd_df$TOP_MODULE <- apply(kd_df[, -c(1)], 1, function(x) {
  colnames(kd_df[, -c(1)])[which.min(x)]
})
# Replace NA with 0 in effect columns
kd_df[is.na(kd_df)] <- 0

network <- read.delim(paste0("Results/complete_analysis/KDA/DEG_SD_metacell/input_files/network/", celltype, "_network.txt"))
network$WEIGHT <- 1
network <- network %>% dplyr::filter(TAIL %in% kd_df$NODE | HEAD %in% kd_df$NODE)

nodes <- data.frame(NODE = unique(c(network$TAIL, network$HEAD)))
nodes <- left_join(nodes, celltype_modfile, by = "NODE")
colnames(nodes) <- gsub(paste0("_", celltype), "", colnames(nodes))

effect_colors <- effect_rgb_colors
kd_df <- kd_df %>%
  mutate(TOP_MODULE_COLOR = effect_colors[TOP_MODULE])

kd_df$colors <- apply(kd_df[, colnames(kd_df) %in% effects], 1, function(row) {
  paste(sprintf("'%s'", effect_colors[names(row)[row > 0]]), collapse = ",") })

kd_df$module_count <- sapply(kd_df$colors, function(x) {
  # Count the number of 'rgb(' in the string
  count <- str_count(x, "rgb\\(")
  # Create a string of '1' repeated based on the count
  paste(rep(1, count), collapse = ",")
})

nodes$colors <- apply(nodes[, colnames(nodes) %in% effects], 1, function(row) {
  paste(sprintf("'%s'", effect_colors[names(row)[row > 0]]), collapse = ",") })
nodes$module_count <- sapply(nodes$colors, function(x) {
  # Count the number of 'rgb(' in the string
  count <- str_count(x, "rgb\\(")
  # Create a string of '1' repeated based on the count
  paste(rep(1, count), collapse = ",")
})
nodes$module_count <- ifelse(nodes$module_count == "", 0, nodes$module_count)

nodes$url <- paste0("https://quickchart.io/chart?c={type:'pie',data:{labels:[],datasets:[{data:[", nodes$module_count, "],backgroundColor:[", nodes$colors, "],},],},options:{legend:{display:false},plugins:{datalabels:{display:false}},elements:{arc:{backgroundColor:'transparent',borderColor:'transparent'},},},}")

nodes[nodes$module_count == 0, "url"] <-"https://quickchart.io/chart?c={type:'pie',data:{labels:[],datasets:[{data:[1],backgroundColor:['rgb(255,255,255)'],},],},options:{legend:{display:false},plugins:{datalabels:{display:false}},elements:{arc:{backgroundColor:'transparent',borderColor:'transparent'},},},}"

kd_df$url <- paste0("https://quickchart.io/chart?c={type:'pie',data:{labels:[],datasets:[{data:[", kd_df$module_count, "],backgroundColor:[", kd_df$colors, "],},],},options:{legend:{display:false},plugins:{datalabels:{display:false}},elements:{arc:{backgroundColor:'transparent',borderColor:'transparent'},},},}")

nodes <- nodes %>% 
  mutate(LABEL = NODE,
         TOP_KD = ifelse(NODE %in% kd_df$NODE, "YES", "NO"),
         SHAPE = ifelse(TOP_KD %in% c("YES"), "Diamond", "Ellipse"),
         LABEL_SIZE = plyr::mapvalues(TOP_KD, 
                                      from = c("NO", "YES"),
                                      to = c("40", "60")),
         WIDTH = plyr::mapvalues(TOP_KD, 
                                 from = c("NO", "YES"),
                                 to = c(100, 200)),
         HEIGHT = WIDTH)

nodes <- nodes %>% 
  left_join(kd_df[, colnames(kd_df) %in% c("NODE", "TOP_MODULE", "TOP_MODULE_COLOR", "module_count", "url")], by = "NODE", suffix = c("", "_new")) %>% 
  mutate(url = coalesce(url_new, url),
         module_count = coalesce(module_count_new, module_count)) %>%             # Replace with new values if available
  dplyr::select(-url_new, -module_count_new)                     

nodes <- nodes %>% dplyr::select(NODE, TOP_KD, LABEL, SHAPE, LABEL_SIZE, WIDTH, HEIGHT, TOP_MODULE_COLOR, url) 

write.table(network, paste0("Results/complete_analysis/KDA/DEG_SD_metacell/sexeffect_separate/", celltype, "/", celltype, "_network_edited.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
write.table(nodes, paste0("Results/complete_analysis/KDA/DEG_SD_metacell/sexeffect_separate/", celltype, "/", celltype, "_nodes_edited.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
