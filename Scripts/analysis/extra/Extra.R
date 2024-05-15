##Extra plots 

##cell type proportion tests

library("scProportionTest")
prop_test <- sc_utils(seurat)

prop_test <- permutation_test(
  prop_test, cluster_identity = "cell.type",
  sample_1 = "XYYM", sample_2 = "XYM",
  sample_identity = "Genotype"
)
plot <- permutation_plot(prop_test)

pdf("Plots/DEG/cell_proportion/XYYM_XYM.pdf",width  = 12, height = 5)
print(plot)
dev.off()

prop_test <- permutation_test(
  prop_test, cluster_identity = "cell.type",
  sample_1 = "XYYF", sample_2 = "XYF",
  sample_identity = "Genotype"
)

plot <- permutation_plot(prop_test)

pdf("Plots/DEG/cell_proportion/XYYF_XYF.pdf",width  = 12, height = 5)
print(plot)
dev.off()

prop_test <- permutation_test(
  prop_test, cluster_identity = "cell.type",
  sample_1 = "XXYM", sample_2 = "XXM",
  sample_identity = "Genotype"
)

plot <- permutation_plot(prop_test)

pdf("Plots/DEG/cell_proportion/XXYM_XXM.pdf",width  = 12, height = 5)
print(plot)
dev.off()

prop_test <- permutation_test(
  prop_test, cluster_identity = "cell.type",
  sample_1 = "XXYF", sample_2 = "XXF",
  sample_identity = "Genotype"
)

plot <- permutation_plot(prop_test)

pdf("Plots/DEG/cell_proportion/XYYF_XXF.pdf",width  = 12, height = 5)
print(plot)
dev.off()

prop_test <- permutation_test(
  prop_test, cluster_identity = "cell.type",
  sample_1 = "XXYM", sample_2 = "XYM",
  sample_identity = "Genotype"
)

plot <- permutation_plot(prop_test)

pdf("Plots/DEG/cell_proportion/XXYM_XYM.pdf",width  = 12, height = 5)
print(plot)
dev.off()

prop_test <- permutation_test(
  prop_test, cluster_identity = "cell.type",
  sample_1 = "XXYF", sample_2 = "XYF",
  sample_identity = "Genotype"
)

plot <- permutation_plot(prop_test)

pdf("Plots/DEG/cell_proportion/XXYF_XYF.pdf",width  = 12, height = 5)
print(plot)
dev.off()

##############################################################################

##heatmap of quadratic genes - finish later 

yquad_genes <- DEG_df %>% dplyr::filter(comparison == "Ydose^2" & adj.P.Val < 0.05) %>% dplyr::select(gene) %>% unlist(use.names = FALSE)
ylin_genes <- DEG_df %>% dplyr::filter(comparison == "Ydose" & adj.P.Val < 0.05) %>% dplyr::select(gene) %>% unlist(use.names = FALSE)

table(yquad_genes %in% ylin_genes)
quad_genes <- yquad_genes[!yquad_genes %in% ylin_genes]

seurat <- ScaleData(seurat, features = rownames(seurat))

DoHeatmap(seurat, features = quad_genes, group.by = Genotype)

##############################################################################

##p-value density plots

##FCG
test_df <- DEG_df_annotated %>% dplyr::filter(model %in% c("Normative", "FCG") & adj.P.Val < 0.05) %>% dplyr::mutate(comparison = factor(comparison, levels = c("XXF_vs_XYM", "Ovary_vs_Testis", "XX_vs_XY", "Ovary:SexChromXX")))

densityplot(~adj.P.Val,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~P.Value | cell,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~logFC,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

densityplot(~logFC | cell,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

##Aneuploidy
test_df <- DEG_df_annotated %>% dplyr::filter(model %in% c("Normative", "Full_aneuploidy") & adj.P.Val < 0.05) %>% dplyr::mutate(comparison = factor(comparison, levels = c("XXF_vs_XYM", "Ovary_vs_Testis", "XX_vs_XY", "Aneuploidy", "Ovary:SexChromXX", "Ovary:Aneuploidy", "SexChromXX:Aneuploidy")))

densityplot(~adj.P.Val,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~adj.P.Val | cell,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~logFC,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

densityplot(~logFC | cell,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

##Full
test_df <- DEG_df_annotated %>% dplyr::filter(model %in% c("Normative", "Full_continuous") & adj.P.Val < 0.05) %>% dplyr::mutate(comparison = factor(comparison, levels = c("XXF_vs_XYM", "Ovary_vs_Testis", "Xdose", "Ydose", "Ovary:Xdose", "Ovary:Ydose", "Xdose:Ydose")))

densityplot(~adj.P.Val,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~adj.P.Val | cell,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~logFC,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

densityplot(~logFC | cell,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

##Quadratic
test_df <- DEG_df_annotated %>% dplyr::filter(model %in% c("Normative", "Full_quadratic") & adj.P.Val < 0.05) %>% dplyr::mutate(comparison = factor(comparison, levels = c("XXF_vs_XYM", "Ovary_vs_Testis", "Xdose", "Ydose", "Ydose^2")))

densityplot(~adj.P.Val,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~adj.P.Val | cell,data=test_df,
            groups=comparison,
            xlab="adj.P.Val",
            main="Adjusted P-Value Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE)

densityplot(~logFC,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

densityplot(~logFC | cell,data=test_df,
            groups=comparison,
            xlab="logFC",
            main="logFC Distribution by Term and Cell Type",
            plot.points=FALSE,
            auto.key=TRUE,
            xlim = c(-2,2))

##############################################################################

library(Seurat)
library(irlba)
library(RSpectra)
install.packages("irlba", type = "source")
install.packages("RSpectra")
source("~/project-xyang123/cerebellum_sct/scUtils/Basic_pipeline/scFunctions.R")

seurat <- readRDS("/u/home/j/jshin/project-xyang123/cerebellum_sct/Data/Derived/PreFilter/seurat.RDS")

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat, features = rownames(seurat))

seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.3)

seurat.subsampled <- seurat[, sample(colnames(seurat), size = 3000, replace=F)]
seurat.subsampled <- RunUMAP(seurat.subsampled, dims = 1:15)

DimPlot(seurat.subsampled, reduction = "umap")
