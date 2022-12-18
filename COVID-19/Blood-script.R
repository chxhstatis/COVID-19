library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/control/predata")
folders <- list.files('./', pattern='a-')
scList <- lapply(folders, function(folder){
  CreateSeuratObject(counts=Read10X(folder),
                     project="control1",
                     min.cells=3, min.features=200)
})

folders <- list.files('./', pattern='b-')
scList2 <- lapply(folders, function(folder){
  CreateSeuratObject(counts=Read10X(folder),
                     project="control2",
                     min.cells=3, min.features=200)
})
rm(folders)
setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")
sample22 <- Read10X(data.dir = "b-sample22")
sample22 <- CreateSeuratObject(counts=sample22, project="patient21",
                               min.cells=3, min.features=200)
sample23 <- Read10X(data.dir = "b-sample23")
sample23 <- CreateSeuratObject(counts=sample23, project="patient22",
                               min.cells=3, min.features=200)
sample24 <- Read10X(data.dir = "b-sample24")
sample24 <- CreateSeuratObject(counts=sample24, project="patient23",
                               min.cells=3, min.features=200)
sample25 <- Read10X(data.dir = "b-sample25")
sample25 <- CreateSeuratObject(counts=sample25, project="patient24",
                               min.cells=3, min.features=200)
sample26 <- Read10X(data.dir = "b-sample26")
sample26 <- CreateSeuratObject(counts=sample26, project="patient25",
                               min.cells=3, min.features=200)
sample27 <- Read10X(data.dir = "b-sample27")
sample27 <- CreateSeuratObject(counts=sample27, project="patient26",
                               min.cells=3, min.features=200)
sample28 <- Read10X(data.dir = "b-sample28")
sample28 <- CreateSeuratObject(counts=sample28, project="patient27",
                               min.cells=3, min.features=200)
sample29 <- Read10X(data.dir = "b-sample29")
sample29 <- CreateSeuratObject(counts=sample29, project="patient28",
                               min.cells=3, min.features=200)
sample30 <- Read10X(data.dir = "b-sample30")
sample30 <- CreateSeuratObject(counts=sample30, project="patient29",
                               min.cells=3, min.features=200)
sample31 <- Read10X(data.dir = "b-sample31")
sample31 <- CreateSeuratObject(counts=sample31, project="patient30",
                               min.cells=3, min.features=200)
sample32 <- Read10X(data.dir = "b-sample32")
sample32 <- CreateSeuratObject(counts=sample32, project="patient31",
                               min.cells=3, min.features=200)
sample33 <- Read10X(data.dir = "b-sample33")
sample33 <- CreateSeuratObject(counts=sample33, project="patient32",
                               min.cells=3, min.features=200)
sample34 <- Read10X(data.dir = "b-sample34")
sample34 <- CreateSeuratObject(counts=sample34, project="patient33",
                               min.cells=3, min.features=200)
sample35 <- Read10X(data.dir = "b-sample35")
sample35 <- CreateSeuratObject(counts=sample35, project="patient34",
                               min.cells=3, min.features=200)
sample36 <- Read10X(data.dir = "b-sample36")
sample36 <- CreateSeuratObject(counts=sample36, project="patient35",
                               min.cells=3, min.features=200)
sample37 <- Read10X(data.dir = "b-sample37")
sample37 <- CreateSeuratObject(counts=sample37, project="patient36",
                               min.cells=3, min.features=200)
sample38 <- Read10X(data.dir = "b-sample38")
sample38 <- CreateSeuratObject(counts=sample38, project="patient37",
                               min.cells=3, min.features=200)
sample39 <- Read10X(data.dir = "b-sample39")
sample39 <- CreateSeuratObject(counts=sample39, project="patient38",
                               min.cells=3, min.features=200)
sample40 <- Read10X(data.dir = "b-sample40")
sample40 <- CreateSeuratObject(counts=sample40, project="patient39",
                               min.cells=3, min.features=200)
sample41 <- Read10X(data.dir = "b-sample41")
sample41 <- CreateSeuratObject(counts=sample41, project="patient40",
                               min.cells=3, min.features=200)
sample42 <- Read10X(data.dir = "b-sample42")
sample42 <- CreateSeuratObject(counts=sample42, project="patient41",
                               min.cells=3, min.features=200)
sample43 <- Read10X(data.dir = "b-sample43")
sample43 <- CreateSeuratObject(counts=sample43, project="patient42",
                               min.cells=3, min.features=200)
sample44 <- Read10X(data.dir = "b-sample44")
sample44 <- CreateSeuratObject(counts=sample44, project="patient43",
                               min.cells=3, min.features=200)
sample45 <- Read10X(data.dir = "b-sample45")
sample45 <- CreateSeuratObject(counts=sample45, project="patient44",
                               min.cells=3, min.features=200)
sample46 <- Read10X(data.dir = "b-sample46")
sample46 <- CreateSeuratObject(counts=sample46, project="patient45",
                               min.cells=3, min.features=200)
sample47 <- Read10X(data.dir = "b-sample47")
sample47 <- CreateSeuratObject(counts=sample47, project="patient46",
                               min.cells=3, min.features=200)
sample48 <- Read10X(data.dir = "b-sample48")
sample48 <- CreateSeuratObject(counts=sample48, project="patient47",
                               min.cells=3, min.features=200)
sample49 <- Read10X(data.dir = "b-sample49")
sample49 <- CreateSeuratObject(counts=sample49, project="patient48",
                               min.cells=3, min.features=200)
sample50 <- Read10X(data.dir = "b-sample50")
sample50 <- CreateSeuratObject(counts=sample50, project="patient49",
                               min.cells=3, min.features=200)
sample52 <- Read10X(data.dir = "b-sample52")
sample52 <- CreateSeuratObject(counts=sample52, project="patient51",
                               min.cells=3, min.features=200)
sample53 <- Read10X(data.dir = "b-sample53")
sample53 <- CreateSeuratObject(counts=sample53, project="patient52",
                               min.cells=3, min.features=200)
sample54 <- Read10X(data.dir = "b-sample54")
sample54 <- CreateSeuratObject(counts=sample54, project="patient53",
                               min.cells=3, min.features=200)
sample55 <- Read10X(data.dir = "b-sample55")
sample55 <- CreateSeuratObject(counts=sample55, project="patient54",
                               min.cells=3, min.features=200)
sample56 <- Read10X(data.dir = "b-sample56")
sample56 <- CreateSeuratObject(counts=sample56, project="patient55",
                               min.cells=3, min.features=200)
sample57 <- Read10X(data.dir = "b-sample57")
sample57 <- CreateSeuratObject(counts=sample57, project="patient56",
                               min.cells=3, min.features=200)
sample60 <- Read10X(data.dir = "b-sample60")
sample60 <- CreateSeuratObject(counts=sample60, project="patient59",
                               min.cells=3, min.features=200)
sample62 <- Read10X(data.dir = "b-sample62")
sample62 <- CreateSeuratObject(counts=sample62, project="patient61",
                               min.cells=3, min.features=200)
sample63 <- Read10X(data.dir = "b-sample63")
sample63 <- CreateSeuratObject(counts=sample63, project="patient62",
                               min.cells=3, min.features=200)
sample64 <- Read10X(data.dir = "b-sample64")
sample64 <- CreateSeuratObject(counts=sample64, project="patient63",
                               min.cells=3, min.features=200)
sample65 <- Read10X(data.dir = "b-sample65")
sample65 <- CreateSeuratObject(counts=sample65, project="patient64",
                               min.cells=3, min.features=200)
Blood <- merge(sample22, c(scList,scList2,sample23,sample24,sample25,
                       sample26,sample27,sample28,sample29,sample30,
                       sample31,sample32,sample33,sample34,sample35,
                       sample36,sample37,sample38,sample39,sample40,
                       sample41,sample42,sample43,sample44,sample45,
                       sample46,sample47,sample48,sample49,sample50,
                       sample52,sample53,sample54,sample55,
                       sample56,sample57,sample60,
                       sample62,sample63,sample64,sample65
                       ))
save(Blood, file = "Blood.Rda")
rm(list = ls())
#-------------------------------------------------------------------------
load("Blood.Rda")
#计算线粒体QC指标
Blood[['percent.mt']] <- PercentageFeatureSet(Blood, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM",
              "HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Blood@assays$RNA)) 
HB.genes <- rownames(Blood@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Blood[["percent.HB"]] <- PercentageFeatureSet(Blood, features=HB.genes)
#核糖体
rb.gene <- rownames(Blood)[grep("^RP[SL]",rownames(Blood))]
Blood <- Blood[!rownames(Blood) %in% rb.gene,]
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
Blood = Blood %>% RunHarmony(group.by.vars="order", plot_convergence = TRUE)
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
ggsave("orig.png", p1+p3, width=40, height=25)
ggsave("umap_tsne.png", p2+p4, width=40, height=25)
save(Blood, file = "Blood2.Rda")
rm(list = ls())
#------------------------------------------------------------------------
load("Blood2.Rda")
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
FeaturePlot(Blood, features = c("LTF","LCN2","DEFA3",
                                "CD3D","CD3E","CD3G","MKI67","NKG7","PRF1","GNLY","KLRB1",
                                "TUBB1","PPBP",
                                "IGHA1","JCHAIN",
                                "TCF4","IRF8",
                                "FCGR3B","S100A8","S100A9","NAMPT",
                                "CD14","VCAN","FCN1","FCGR3A",
                                "CLC","GATA2","FCER1A",
                                "HLA-DPA1","HLA-DPB1","CD1C",
                                "MS4A1","CD79A","CD79B"),                     
            label = T, order = T)
VlnPlot(Blood, features=c("LTF","LCN2","DEFA3",
                          "CD3D","CD3E","CD3G","MKI67","NKG7","PRF1","GNLY","KLRB1",
                          "TUBB1","PPBP",
                          "IGHA1","JCHAIN",
                          "TCF4","IRF8",
                          "FCGR3B","S100A8","S100A9","NAMPT",
                          "CD14","VCAN","FCN1","FCGR3A",
                          "CLC","GATA2","FCER1A",
                          "HLA-DPA1","HLA-DPB1","CD1C",
                          "MS4A1","CD79A","CD79B"),
        pt.size=0)
DotPlot(Blood, features = c("LTF","LCN2","DEFA3",
                            "CD3D","CD3E","CD3G","MKI67","NKG7","PRF1","GNLY","KLRB1",
                            "TUBB1","PPBP",
                            "IGHA1","JCHAIN",
                            "TCF4","IRF8",
                            "FCGR3B","S100A8","S100A9","NAMPT",
                            "CD14","VCAN","FCN1","FCGR3A",
                            "CLC","GATA2","FCER1A",
                            "HLA-DPA1","HLA-DPB1","CD1C",
                            "MS4A1","CD79A","CD79B"))
#自定义大群注释
new.cluster.ids <- c("Neutrophils","Neutrophils","T cells","Neutrophils","Neutrophils",#0-4
                     "T cells","NK cells","CD14 Monocytes","CD14 Monocytes","T cells",      #5-9  
                     "Platelets","B cells","Neutrophils","CD14 Monocytes","NK cells",  #10-14     
                     "CD14 Monocytes","Neutrophils","Neutrophils","FCGR3A+ Monocytes","Platelets",#15-19
                     "Transferrin-associated neutrophils","T cells","Proliferating NK/T cells","B cells","DCs",  #20-24
                     "Plasma cells","T cells","Eosinophils","pDCs","DCs",#25-29
                     "pDCs","NK cells","Plasma cells")#30-32
names(new.cluster.ids) <- levels(Blood)
Blood <- RenameIdents(Blood, new.cluster.ids)
Blood[["celltype"]] <- Idents(Blood)
p2 <- DimPlot(Blood, reduction = "tsne", label = TRUE, pt.size = 0.5)
p1 <- DimPlot(Blood, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 + p2
save(Blood, file = "Blood3.Rda")

