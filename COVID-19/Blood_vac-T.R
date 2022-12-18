library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)

load("Bloodvac3.Rda")

#取T细胞
T_cell <- subset(Blood, idents = c("T cells"))
#标准化
T_cell <- NormalizeData(T_cell, normalization.method = "LogNormalize",
                        scale.factor = 10000)
#识别高度可变的特征（特征选择）
T_cell <- FindVariableFeatures(T_cell, selection.method = "vst",
                               nfeatures = 2000)
#缩放数据(all.genes)
all.genes <- rownames(T_cell)
T_cell <- ScaleData(T_cell, features = all.genes)
#PCA
T_cell <- RunPCA(T_cell, npcs = 50, 
                 features = VariableFeatures(object = T_cell),
                 verbose = F)
DimPlot(T_cell, reduction = "pca")+ NoLegend()
ElbowPlot(T_cell)
#harmony矫正
T_cell = T_cell %>% RunHarmony(group.by.vars="day", plot_convergence = TRUE)
DimPlot(T_cell, reduction = "harmony")+ NoLegend()
#细胞聚类
T_cell <- T_cell %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
T_cell <- T_cell %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
#降维可视化
p1 <- DimPlot(T_cell, reduction = "tsne", group.by = "orig.ident",
              pt.size = 0.5)
p2 <- DimPlot(T_cell, reduction = "tsne", pt.size = 0.5, label=T)
p3 <- DimPlot(T_cell, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.5)
p4 <- DimPlot(T_cell, reduction = "umap", pt.size = 0.5, label=T)
p1+p3
p2+p4
save(T_cell, file = "T_cell.Rda")
rm(list = ls())

#-------------------------------------------------------------------
load("T_cell.Rda")
T_cell.markers <- FindAllMarkers(T_cell, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- T_cell.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top10_markers <- T_cell.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
save(T_cell.markers, top10_markers,
     file = "T_cell_marker.Rda")
VlnPlot(T_cell, features=c("SLC4A10","CD3D","CD3E",
                           "ISG15","STAT1","LY6E","OAS1",
                           "FCER1G","NAMPT","S100A8","IFITM1","IFITM2","IFITM3",
                           "TRDV2","TRGV9",
                           "CD8A","CD4",
                           "SELL","CCR7","TCF7",
                           "LAG3","PDCD1","HAVCR2","TBX21",
                           "CXCR4","CXCR3","CD44","GZMK","GZMB",
                           "PRF1","CX3CR1","CCL5",
                           "IL2RA","FOXP3",
                           "CTLA4","CCR4","CD69","IL7R","LTB"),
        pt.size=0)
FeaturePlot(T_cell, features=c("SLC4A10","CD3D","CD3E",
                               "ISG15","STAT1","LY6E","OAS1",
                               "FCER1G","NAMPT","S100A8","IFITM1","IFITM2","IFITM3",
                               "TRDV2","TRGV9",
                               "CD8A","CD4",
                               "SELL","CCR7","TCF7",
                               "LAG3","PDCD1","HAVCR2","TBX21",
                               "CXCR4","CXCR3","CD44","GZMK","GZMB",
                               "PRF1","CX3CR1","CCL5",
                               "IL2RA","FOXP3",
                               "CTLA4","CCR4","CD69","IL7R","LTB"),
            label=T,order=T)
DotPlot(T_cell, features=c("SLC4A10","CD3D","CD3E",
                           "ISG15","STAT1","LY6E","OAS1",
                           "FCER1G","NAMPT","S100A8","IFITM1","IFITM2","IFITM3",
                           "TRDV2","TRGV9",
                           "CD8A","CD4",
                           "SELL","CCR7","TCF7",
                           "LAG3","PDCD1","HAVCR2","TBX21",
                           "CXCR4","CXCR3","CD44","GZMK","GZMB",
                           "PRF1","CX3CR1","CCL5",
                           "IL2RA","FOXP3",
                           "CTLA4","CCR4","CD69","IL7R","LTB"))

new.cluster.ids <- c("CD8Teff","CD4Tn","Resting CD4T","Activated CD4T","CD4Tn",
                     "CD4Tm","Resting CD4T","CD8Tem","CD8Tn","CD8Teff",           
                     "Resting CD4T")
names(new.cluster.ids) <- levels(T_cell)
T_cell <- RenameIdents(T_cell, new.cluster.ids)
T_cell[["celltype"]] <- Idents(T_cell)
save(T_cell, file = "T_cell2.Rda")
p2 <- DimPlot(T_cell, reduction = "tsne", label = TRUE, pt.size = 0.5)
p1 <- DimPlot(T_cell, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 + p2