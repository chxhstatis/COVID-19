library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")

##innate-like Tcells
load('T_cell2.Rda')
T15 <- subset(T_cell,idents=c('innate-like Tcells'))

#标准化
T15 <- NormalizeData(T15, normalization.method = "LogNormalize",
                     scale.factor = 10000)
#识别高度可变的特征（特征选择）
T15 <- FindVariableFeatures(T15, selection.method = "vst",
                            nfeatures = 2000)
#缩放数据(all.genes)
all.genes <- rownames(T15)
T15 <- ScaleData(T15, features = all.genes)
#PCA
T15 <- RunPCA(T15, npcs = 50, 
              features = VariableFeatures(object = T15),
              verbose = F)
DimPlot(T15, reduction = "pca")+ NoLegend()
ElbowPlot(T15)
#harmony矫正
T15 = T15 %>% RunHarmony(group.by.vars="order", plot_convergence = TRUE)
DimPlot(T15, reduction = "harmony")+ NoLegend()
#细胞聚类
T15 <- T15 %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2) %>%
  identity()
T15 <- T15 %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2) %>%
  identity()
#降维可视化
p1 <- DimPlot(T15, reduction = "tsne", group.by = "order",
              pt.size = 0.5)
p2 <- DimPlot(T15, reduction = "tsne", pt.size = 0.5, label=T)
p3 <- DimPlot(T15, reduction = "umap", group.by = "order",
              pt.size = 0.5)
p4 <- DimPlot(T15, reduction = "umap", pt.size = 0.5, label=T)
p1+p3
p2+p4

new.cluster.ids <- c("C0","C1")
names(new.cluster.ids) <- levels(T15)
T15 <- RenameIdents(T15, new.cluster.ids)
T15[["celltype"]] <- Idents(T_cell)
save(T15, file = "T15.Rda")

T15.markers <- FindAllMarkers(T15, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
top10_markers <- T15.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

save(T_cell.markers, top10_markers,
     file = "T15_marker.Rda")