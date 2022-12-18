library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/control/predata")
folders <- list.files('./', pattern='b-')
scList <- lapply(folders, function(folder){
  CreateSeuratObject(counts=Read10X(folder),
                     project="control",
                     min.cells=3, min.features=200)
})
rm(folders)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/vaccine")
control1 <- readRDS("control7.rds")
control1 <- control1$exon
control1 <- CreateSeuratObject(counts=control1, project="control",
                               min.cells=3, min.features=200)
control2 <- readRDS("control8.rds")
control2 <- control2$exon
control2 <- CreateSeuratObject(counts=control2, project="control",
                               min.cells=3, min.features=200)
control3 <- readRDS("control9.rds")
control3 <- control3$exon
control3 <- CreateSeuratObject(counts=control3, project="control",
                               min.cells=3, min.features=200)
control4 <- readRDS("control10.rds")
control4 <- control4$exon
control4 <- CreateSeuratObject(counts=control4, project="control",
                               min.cells=3, min.features=200)
control5 <- readRDS("control11.rds")
control5 <- control5$exon
control5 <- CreateSeuratObject(counts=control5, project="control",
                               min.cells=3, min.features=200)
control6 <- readRDS("control12.rds")
control6 <- control6$exon
control6 <- CreateSeuratObject(counts=control6, project="control",
                               min.cells=3, min.features=200)
vaccine1 <- Read10X(data.dir = "BBday1")
vaccine1 <- CreateSeuratObject(counts=vaccine1, project="day1",
                               min.cells=3, min.features=200)
vaccine2 <- Read10X(data.dir = "BBday7")
vaccine2 <- CreateSeuratObject(counts=vaccine2, project="day7",
                               min.cells=3, min.features=200)
vaccine3 <- Read10X(data.dir = "BB2day1")
vaccine3 <- CreateSeuratObject(counts=vaccine3, project="day1",
                               min.cells=3, min.features=200)
vaccine4 <- Read10X(data.dir = "BB2day7")
vaccine4 <- CreateSeuratObject(counts=vaccine4, project="day7",
                               min.cells=3, min.features=200)
vaccine5 <- Read10X(data.dir = "CCday1")
vaccine5 <- CreateSeuratObject(counts=vaccine5, project="day1",
                               min.cells=3, min.features=200)
vaccine6 <- Read10X(data.dir = "CCday7")
vaccine6 <- CreateSeuratObject(counts=vaccine6, project="day7",
                               min.cells=3, min.features=200)
vaccine7 <- Read10X(data.dir = "CC2day1")
vaccine7 <- CreateSeuratObject(counts=vaccine7, project="day1",
                               min.cells=3, min.features=200)
vaccine8 <- Read10X(data.dir = "CC2day7")
vaccine8 <- CreateSeuratObject(counts=vaccine8, project="day7",
                               min.cells=3, min.features=200)
vaccine9 <- Read10X(data.dir = "CBday1")
vaccine9 <- CreateSeuratObject(counts=vaccine9, project="day1",
                               min.cells=3, min.features=200)
vaccine10 <- Read10X(data.dir = "CBday7")
vaccine10 <- CreateSeuratObject(counts=vaccine10, project="day7",
                               min.cells=3, min.features=200)
vaccine11 <- Read10X(data.dir = "CB2day1")
vaccine11 <- CreateSeuratObject(counts=vaccine11, project="day1",
                               min.cells=3, min.features=200)
vaccine12 <- Read10X(data.dir = "CB2day7")
vaccine12 <- CreateSeuratObject(counts=vaccine12, project="day7",
                               min.cells=3, min.features=200)

Blood <- merge(control1, c(control2,control3,control4,control5,control6,
                           vaccine1,vaccine2,vaccine3,vaccine4,
                           vaccine5,vaccine6,vaccine7,vaccine8,vaccine9,
                           vaccine10,vaccine11,vaccine12))
save(Blood, file = "Bloodvac.Rda")
rm(list = ls())
#-------------------------------------------------------------------------
load("Bloodvac.Rda")
#核糖体
rb.gene <- rownames(Blood)[grep("^RP[SL]",rownames(Blood))]
Blood <- Blood[!rownames(Blood) %in% rb.gene,]
#计算线粒体QC指标
Blood[['percent.mt']] <- PercentageFeatureSet(Blood, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM",
              "HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Blood@assays$RNA)) 
HB.genes <- rownames(Blood@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Blood[["percent.HB"]] <- PercentageFeatureSet(Blood, features=HB.genes)
#可视化QC指标
VlnPlot(Blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                            "percent.HB"), ncol = 4, pt.size=0)
plot1 <- FeatureScatter(Blood, feature1 = "nCount_RNA",
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(Blood, feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
plot1+plot2
rm(plot1, plot2, HB_m, HB.genes, rb.gene)
#过滤
Blood <- subset(Blood, subset = nFeature_RNA>200 & nFeature_RNA<7500 &
                  percent.mt<15 & percent.HB<5)
#标准化
Blood <- NormalizeData(Blood, normalization.method = "LogNormalize",
                       scale.factor = 10000)
#识别高度可变的特征（特征选择）
Blood <- FindVariableFeatures(Blood, selection.method = "vst",
                              nfeatures = 2000)
#缩放数据(all.genes)
all.genes <- rownames(Blood)
Blood <- ScaleData(Blood, features = all.genes)
#PCA
Blood <- RunPCA(Blood, npcs = 50, 
                features = VariableFeatures(object = Blood),
                verbose = F)
DimPlot(Blood, reduction = "pca")+ NoLegend()
ElbowPlot(Blood)
#harmony矫正
Blood = Blood %>% RunHarmony(group.by.vars="day", plot_convergence = TRUE)
DimPlot(Blood, reduction = "harmony")+ NoLegend()
#细胞聚类
Blood <- Blood %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.0) %>%
  identity()
Blood <- Blood %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.0) %>%
  identity()
#降维可视化
p1 <- DimPlot(Blood, reduction = "tsne", group.by = "orig.ident",
              pt.size = 0.5)
p2 <- DimPlot(Blood, reduction = "tsne", pt.size = 0.5, label=T)
p3 <- DimPlot(Blood, reduction = "umap", group.by = "orig.ident",
              pt.size = 0.5)
p4 <- DimPlot(Blood, reduction = "umap", pt.size = 0.5, label=T)
save(Blood, file = "Bloodvac2.Rda")
rm(list = ls())
#------------------------------------------------------------------------
load("Bloodvac=2.Rda")
#寻找cluster差异表达的特征
Blood.markers <- FindAllMarkers(Blood, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- Blood.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top10_markers <- Blood.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
top5_markers
DotPlot(Blood, features = unique(top5_markers$gene)) &
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1)) & NoLegend()
save(Blood.markers, top5_markers,
     file = "Blood_marker.Rda")
FeaturePlot(Blood, features = c("CD3D","CD3E","CD3G","MKI67","NKG7","PRF1","GNLY","KLRB1",
                                "TUBB1","PPBP",
                                "IGHA1","JCHAIN",
                                "TCF4","IRF8",
                                "FCGR3B","S100A8","S100A9","NAMPT",
                                "CD14","VCAN","FCN1","FCGR3A",
                                "CLC","GATA2","FCER1A",
                                "HLA-DPA1","HLA-DPB1","CD1C",
                                "MS4A1","CD79A","CD79B"),                     
            label = T, order = T)
VlnPlot(Blood, features=c("CD3D","CD3E","CD3G","MKI67","NKG7","PRF1","GNLY","KLRB1",
                          "TUBB1","PPBP",
                          "IGHA1","JCHAIN",
                          "TCF4","IRF8",
                          "FCGR3B","S100A8","S100A9","NAMPT",
                          "CD14","VCAN","FCN1","FCGR3A",
                          "CLC","GATA2","FCER1A",
                          "HLA-DPA1","HLA-DPB1","CD1C",
                          "MS4A1","CD79A","CD79B"),
        pt.size=0)
DotPlot(Blood, features = c("CD3D","CD3E","CD3G","MKI67","NKG7","PRF1","GNLY","KLRB1",
                            "TUBB1","PPBP",
                            "IGHA1","JCHAIN",
                            "TCF4","IRF8",
                            "FCGR3B","S100A8","S100A9","NAMPT",
                            "CD14","VCAN","FCN1","FCGR3A",
                            "CLC","GATA2","FCER1A",
                            "HLA-DPA1","HLA-DPB1","CD1C",
                            "MS4A1","CD79A","CD79B"))
#自定义大群注释
new.cluster.ids <- c("T cells","T cells","NK cells","CD14 Monocytes","T cells",#0-4
                     "NK cells","B cells","NK cells","T cells","T cells",      #5-9  
                     "T cells","T cells","B cells","NK cells","CD14 Monocytes",  #10-14     
                     "FCGR3A+ Monocytes","B cells","Platelets","CD14 Monocytes","T cells",#15-19
                     "Proliferating NK/T cells","CD14 Monocytes","DCs","B cells","Plasma cells",  #20-24
                     "pDCs","Others","Others")#25-27
names(new.cluster.ids) <- levels(Blood)
Blood <- RenameIdents(Blood, new.cluster.ids)
Blood[["celltype"]] <- Idents(Blood)
p2 <- DimPlot(Blood, reduction = "tsne", label = TRUE, pt.size = 0.5)
p1 <- DimPlot(Blood, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 + p2
save(Blood, file = "Bloodvac3.Rda")
