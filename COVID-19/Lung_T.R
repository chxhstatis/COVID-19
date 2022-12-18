library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Lung/patient/predata")
load("Lung3.Rda")

T_NK <- subset(Lung, idents = c("T/NK"))
#标准化
T_NK <- NormalizeData(T_NK, normalization.method = "LogNormalize",
                      scale.factor = 10000)
#识别高度可变的特征（特征选择）
T_NK <- FindVariableFeatures(T_NK, selection.method = "vst",
                             nfeatures = 2000)
#缩放数据(all.genes)
all.genes <- rownames(T_NK)
T_NK <- ScaleData(T_NK, features = all.genes)
#PCA
T_NK <- RunPCA(T_NK, npcs = 50, 
               features = VariableFeatures(object = T_NK),
               verbose = F)
DimPlot(T_NK, reduction = "pca")+ NoLegend()
ElbowPlot(T_NK)
#harmony矫正
T_NK = T_NK %>% RunHarmony(group.by.vars="batch", plot_convergence = TRUE)
DimPlot(T_NK, reduction = "harmony")+ NoLegend()
#细胞聚类
T_NK <- T_NK %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2) %>%
  identity()
T_NK <- T_NK %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2) %>%
  identity()
#降维可视化
p1 <- DimPlot(T_NK, reduction = "tsne", group.by = "orig.ident",
              pt.size = 0.5)
p2 <- DimPlot(T_NK, reduction = "tsne", pt.size = 0.5, label=T)
p3 <- DimPlot(T_NK, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.5)
p4 <- DimPlot(T_NK, reduction = "umap", pt.size = 0.5, label=T)
p1+p3
p2+p4
save(T_NK, file = "T_NK_2.Rda")
#寻找cluster差异表达的特征
T_NK.markers <- FindAllMarkers(T_NK, only.pos = TRUE,
                               min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- T_NK.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top10_markers <- T_NK.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DotPlot(T_NK, features = unique(top5_markers$gene)) &
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1)) & NoLegend()
save(T_NK.markers, top10_markers,
     file = "T_NK_marker.Rda")

DotPlot(T_NK, features=c("CD3D","CD4","CD8A","CD8B",
                         "FCGR3A","NKTR","NCAM1",
                         "IL7R","SELL","CCR7","TCF7","CD27",
                         "PRF1","CCL5","GZMH","FASLG","CD44","CX3CR1",
                         "CD28","CTLA4",
                         "TRDV2","TRGV9",
                         "IL2RA","FOXP3","ZNF683",
                         "LAG3","PDCD1","HAVCR2","TNFRSF9"))
FeaturePlot(T_NK, features=c("CD3D","CD4","CD8A","CD8B",
                             "FCGR3A","NKTR","NCAM1",
                             "IL7R","SELL","CCR7","KLRG1","CD27",
                             "PRF1","CCL5","GZMH","FASLG","CD44","CX3CR1",
                             "CD28","CTLA4",
                             "TRDV2","TRGV9",
                             "IL2RA","FOXP3","ZNF683",
                             "LAG3","PDCD1","HAVCR2","TNFRSF9"),
            label=T, order=T)
VlnPlot(T_NK, features=c("CD3D","CD4","CD8A","CD8B",
                         "FCGR3A","NKTR","NCAM1",
                         "IL7R","SELL","CCR7","KLRG1","CD27",
                         "PRF1","CCL5","GZMH","FASLG","CD44","CX3CR1",
                         "CD28","CTLA4",
                         "TRDV2","TRGV9",
                         "IL2RA","FOXP3","ZNF683",
                         "LAG3","PDCD1","HAVCR2","TNFRSF9"),
        pt.size=0)
new.cluster.ids <- c("Resting CD4T","CTL","Activated CD4T", "Resting CD8T","NK cells",
                     "CD4Treg","mixed Tcells","innate-like Tcells","ISG Tcells","CTL")
names(new.cluster.ids) <- levels(T_NK)
T_NK <- RenameIdents(T_NK, new.cluster.ids)
p2 <- DimPlot(T_NK, reduction = "tsne", label = TRUE, pt.size = 0.5)
p1 <- DimPlot(T_NK, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 + p2
T_NK[["celltype"]] <- Idents(T_NK)
save(T_NK,file="T_NK_3.Rda")