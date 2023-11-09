#### Analysis of CAF subtypes using scRNA-Seq data from Peng et al. ####

#Loading Packages
library(Seurat)
library(ggplot2)
library(scales)
library(ggpubr)
library(RColorBrewer)
library(harmony)
library(qusage)

setwd("Working Directory")
fibro <- readRDS("Fibroblasts_peng_et_al.rds")

#Select sample with 100 or more cells
Idents(fibro) <- fibro$Sample
a <- table(Idents(fibro))
fibro <- subset(fibro, idents = names(a[a>=100]))

#Normalization and integration using Harmony
fibro <- SCTransform(fibro,
                     vars.to.regress = "percent.mt",
                     assay = "RNA",
                     return.only.var.genes = T,
                     variable.features.n = 3000)

fibro <- RunPCA(fibro, assay = "SCT", npcs = 50)

fibro <- RunHarmony(fibro,
                    group.by.vars = "Sample",
                    theta = 1,
                    plot_convergence = TRUE)

ElbowPlot(fibro, ndims = 50, reduction = "harmony")

##Evaluate integration
#Compute UMAP before (pca) or after (harmony) integration
fibro <- RunUMAP(object = fibro, reduction = "harmony", dims = 1:20)
fibro <- RunUMAP(object = fibro, reduction = "pca", dims = 1:20, reduction.name = "umap_pca")


plot1 <- DimPlot(fibro, label = T, pt.size = 2, group.by = "Sample") + theme(legend.position = "none")
plot2 <- DimPlot(fibro, label = T, pt.size = 2, reduction = "umap_pca", group.by = "Sample") + theme(legend.position = "none")
plot1+plot2

#Clustering
fibro <- FindNeighbors(object = fibro, reduction = "harmony", dims = 1:20)
fibro <- FindClusters(object = fibro, resolution = 0.5)

#Calculate percentage of ribosomal genes
fibro[["percent.rb"]] <- PercentageFeatureSet(fibro, pattern = "^RP[SL]")

#Plot used for Quality Control
VlnPlot(fibro, stack = T, features = c("LUM", "KRT19", "PTPRC", "PRSS1", "RGS5", "percent.mt", "percent.rb")) +
  geom_boxplot(width = 0.05, fill = "white")

#Selection of clusters to work with
fibro_filtered <-
                subset(fibro, invert = T,
                idents = c(4, #high percent.rb
                           6, #PTPRC/CD45, high percent.mt
                           7, #PRSS1, high percent.mt
                           8, #KRT19, high percent.mt
                           10 #Low LUM, Peri-islet Schwann Cell?
                ))

#Performing again Clustering and dimensional reduction
fibro_filtered <- RunUMAP(object = fibro_filtered, reduction = "harmony", dims = 1:20)
fibro_filtered <- FindNeighbors(object = fibro_filtered, reduction = "harmony", dims = 1:20)
fibro_filtered <- FindClusters(object = fibro_filtered, resolution = 0.25)

#Renaming clusters
fibro_filtered <- 
  RenameIdents(fibro_filtered,
               '0' = "myCAF",
               '1' = "NAF",
               '2' = "Transitional_CAF",
               '3' = "iCAF")


Idents(fibro_filtered) <- factor(as.character(Idents(fibro_filtered)), levels = rev(c("iCAF", "myCAF", "Transitional_CAF", "NAF")))

fibro_filtered$Fibroblast_Clusters <- Idents(fibro_filtered)

#Add score for sialylation pathways
gs <- read.gmt("Sia_gs.gmt")
fibro_filtered <- AddModuleScore(fibro_filtered, features = gs, assay = "SCT")
colnames(fibro_filtered@meta.data)[substr(colnames(fibro_filtered@meta.data), 1, 5) == "Clust"] <- names(gs)

##### PLOTS #####
##Figure
#Violin Plot - CAF subtypes Markers
Idents(fibro_filtered) <- fibro_filtered$Fibroblast_Clusters
VlnPlot(fibro_filtered, stack = T, features = c("LUM", "COL1A1", "COL1A2",
                                       "IGF1", "CCL2", "PLA2G2A", "LMNA",
                                       "COL10A1", "MMP11", "POSTN",
                                       "ACTA2")) +
  geom_boxplot(width = 0.05, fill = "white") +
  theme(legend.position = "none", axis.title.y = element_blank())


#Feature Plot - CAF subtypes Markers
FeaturePlot(fibro_filtered, ncol= 3,
            features = c("FAP", "THY1", "CD34",
                         "MCAM", "NGFR", "ACTA2",
                         "PDPN", "PDGFRA", "VCAM1")) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) &
  theme(legend.position = "none")

##Figures
#Comparing CAF subtypes in Adjacent Normal and Tumor
#Select clusters with at least 100 cells

fibro_filtered$CAF_NT <- paste(Idents(fibro_filtered), fibro_filtered$Condition, sep = "_")
fibro_filtered$CAF_NT <- factor(as.character(fibro_filtered$CAF_NT),
          levels = rev(c("iCAF_T", "myCAF_T", "Transitional_CAF_T", "NAF_T", "NAF_N")))
Idents(fibro_filtered) <- fibro_filtered$CAF_NT
a <- table(Idents(fibro_filtered))
fibro_NT <- subset(fibro_filtered, idents = names(a[a>=100]))

#Violin Plots with pairwise comparisons
data <- fibro_NT@meta.data[,c("Condition", "Fibroblast_Clusters", "CAF_NT", names(gs))]
data2 <- t(as.matrix(fibro_NT@assays$SCT@data[c("ST3GAL4", "NANS", "CMAS", "ST6GALNAC6"),]))
data <- cbind(data, data2)
rm("data2")

ggplot(data, aes(x = CAF_NT, y = ST3GAL4)) +
  geom_violin(aes(fill = CAF_NT), scale = "width") +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_pwc(method = "wilcox_test", label = "p.adj.signif", p.adjust.method = "holm",
           hide.ns = T, step.increase = 0.125, bracket.nudge.y = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")


##Figure
#Density Plots

umap <- fibro_filtered@reductions$umap@cell.embeddings
umap <- as.data.frame(umap)
umap_normal <- umap[fibro_filtered$Condition == "N",] 
umap_tumor <- umap[fibro_filtered$Condition == "T",] 
umap$Condition <- fibro_filtered$Condition

plot1 <- ggplot(umap_normal, aes(x=UMAP_1, y=UMAP_2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 10) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  xlim(c(-12, 12)) +
  ylim(c(-6, 6)) + theme_classic() + theme(legend.position = "none")


plot2 <- ggplot(umap_tumor, aes(x=UMAP_1, y=UMAP_2)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins = 10) +
  scale_fill_gradient2(high = muted("red"), low = muted("blue")) +
  xlim(c(-12, 12)) +
  ylim(c(-6, 6)) + theme_classic() + theme(legend.position = "none")

plot1 + plot2

