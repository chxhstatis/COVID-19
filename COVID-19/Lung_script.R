library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(RColorBrewer)

#数据读入
setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Lung/control/predata")
a <- list.files(pattern=".csv")
data=read.csv(a[1], sep=",", row.names=1)
for(i in 2:10){
  datat=read.csv(a[i], sep=",", row.names=1)
  data=cbind(data,datat)
  rm(datat)
  print(i)
}
data <- as(as.matrix(data), "dgCMatrix")
c1_10 <- CreateSeuratObject(counts=data, project="control",
                            min.cells=3, min.features=200)
rm(data,a,i)
setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Lung/patient/predata")
sample1 <- Read10X(data.dir="p-sample1")
sample1 <- CreateSeuratObject(counts=sample1, project="patient1",
                              min.cells=3, min.features=200)
sample2 <- Read10X(data.dir="p-sample2")
sample2 <- CreateSeuratObject(counts=sample2, project="patient2",
                              min.cells=3, min.features=200)
sample3 <- Read10X(data.dir="p-sample3")
sample3 <- CreateSeuratObject(counts=sample3, project="patient3",
                              min.cells=3, min.features=200)
sample4 <- Read10X(data.dir="p-sample4")
sample4 <- CreateSeuratObject(counts=sample4, project="patient4",
                              min.cells=3, min.features=200)
sample5 <- Read10X(data.dir="p-sample5")
sample5 <- CreateSeuratObject(counts=sample5, project="patient5",
                              min.cells=3, min.features=200)
sample6 <- Read10X(data.dir="p-sample6")
sample6 <- CreateSeuratObject(counts=sample6, project="patient6Alive",
                              min.cells=3, min.features=200)
sample7 <- Read10X(data.dir="p-sample7")
sample7 <- CreateSeuratObject(counts=sample7, project="patient6Dead",
                              min.cells=3, min.features=200)
sample8 <- Read10X(data.dir="p-sample8")
sample8 <- CreateSeuratObject(counts=sample8, project="patient7",
                              min.cells=3, min.features=200)
sample9 <- Read10X(data.dir="p-sample9")
sample9 <- CreateSeuratObject(counts=sample9, project="patient8",
                              min.cells=3, min.features=200)
sample10 <- Read10X(data.dir="p-sample10")
sample10 <- CreateSeuratObject(counts=sample10, project="patient9",
                              min.cells=3, min.features=200)
sample11 <- Read10X(data.dir="p-sample11")
sample11 <- CreateSeuratObject(counts=sample11, project="patient10",
                              min.cells=3, min.features=200)
sample12 <- Read10X(data.dir="p-sample12")
sample12 <- CreateSeuratObject(counts=sample12, project="patient11",
                              min.cells=3, min.features=200)
sample13 <- Read10X(data.dir="p-sample13")
sample13 <- CreateSeuratObject(counts=sample13, project="patient12",
                              min.cells=3, min.features=200)
sample14 <- Read10X(data.dir="p-sample14")
sample14 <- CreateSeuratObject(counts=sample14, project="patient13",
                              min.cells=3, min.features=200)
sample15 <- Read10X(data.dir="p-sample15")
sample15 <- CreateSeuratObject(counts=sample15, project="patient14",
                              min.cells=3, min.features=200)
sample16 <- Read10X(data.dir="p-sample16")
sample16 <- CreateSeuratObject(counts=sample16, project="patient15",
                              min.cells=3, min.features=200)
sample17 <- Read10X(data.dir="p-sample17")
sample17 <- CreateSeuratObject(counts=sample17, project="patient16",
                              min.cells=3, min.features=200)
sample18 <- Read10X(data.dir="p-sample18")
sample18 <- CreateSeuratObject(counts=sample18, project="patient17",
                              min.cells=3, min.features=200)
sample19 <- Read10X(data.dir="p-sample19")
sample19 <- CreateSeuratObject(counts=sample19, project="patient18",
                              min.cells=3, min.features=200)
sample20 <- Read10X(data.dir="p-sample20")
sample20 <- CreateSeuratObject(counts=sample20, project="patient19",
                              min.cells=3, min.features=200)
sample21 <- Read10X(data.dir="p-sample21")
sample21 <- CreateSeuratObject(counts=sample21, project="patient20",
                              min.cells=3, min.features=200)
Lung <- merge(c1_10, c(sample1,sample2,sample3,sample4,
                       sample5,sample6,sample7,sample8,sample9,
                       sample10,sample11,sample12,sample13,sample14,
                       sample15,sample16,sample17,sample18,sample19,
                       sample20,sample21))
save(Lung, file="Lung.Rda")
rm(list = ls())
#--------------------------------------------------------------------
load("Lung.Rda")
#核糖体
rb.gene <- rownames(Lung)[grep("^RP[SL]",rownames(Lung))]
Lung <- Lung[!rownames(Lung) %in% rb.gene,]
#计算线粒体QC指标
Lung[['percent.mt']] <- PercentageFeatureSet(Lung, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM",
              "HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Lung@assays$RNA)) 
HB.genes <- rownames(Lung@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
Lung[["percent.HB"]] <- PercentageFeatureSet(Lung, features=HB.genes)
#可视化QC指标
VlnPlot(Lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                            "percent.HB"), ncol = 4, pt.size=0)
plot1 <- FeatureScatter(Lung, feature1 = "nCount_RNA",
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(Lung, feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
plot1+plot2
rm(plot1, plot2, HB_m, HB.genes,rb.gene)
#过滤
Lung <- subset(Lung, subset = nFeature_RNA>200 & nFeature_RNA<7500 &
                  percent.mt<15 & percent.HB<5)
#标准化
Lung <- NormalizeData(Lung, normalization.method = "LogNormalize",
                       scale.factor = 10000)
#识别高度可变的特征（特征选择）
Lung <- FindVariableFeatures(Lung, selection.method = "vst",
                              nfeatures = 2000)
#缩放数据(all.genes)
all.genes <- rownames(Lung)
Lung <- ScaleData(Lung, features = all.genes)
#PCA
Lung <- RunPCA(Lung, npcs = 50, 
                features = VariableFeatures(object = Lung),
                verbose = F)
DimPlot(Lung, reduction = "pca")+ NoLegend()
ElbowPlot(Lung)
#harmony矫正
Lung = Lung %>% RunHarmony(group.by.vars="batch", plot_convergence = TRUE)
DimPlot(Lung, reduction = "harmony")+ NoLegend()
#细胞聚类
Lung <- Lung %>%
  RunTSNE(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.0) %>%
  identity()
Lung <- Lung %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1.0) %>%
  identity()
#降维可视化
p1 <- DimPlot(Lung, reduction = "tsne", group.by = "orig.ident",
              pt.size = 0.5)
p2 <- DimPlot(Lung, reduction = "tsne", pt.size = 0.5, label=T)
p3 <- DimPlot(Lung, reduction = "umap", group.by = "batch",
              pt.size = 0.5)
p4 <- DimPlot(Lung, reduction = "umap", pt.size = 0.5, label=T)
ggsave("orig.png", p1+p3, width=25, height=15)
ggsave("umap_tsne.png", p2+p4, width=25, height=15)
save(Lung, file = "Lung2.Rda")
rm(list = ls())
#------------------------------------------------------------------------
load("Lung2.Rda")
#寻找cluster差异表达的特征
Lung.markers <- FindAllMarkers(Lung, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25)
top5_markers <- Lung.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top10_markers <- Lung.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
top5_markers
DotPlot(Lung, features = unique(top5_markers$gene)) &
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 8, hjust = 1)) & NoLegend()
ggsave("top5.png", c, width=45, height=20)
save(Lung.markers, top5_markers,
     file = "Lung_marker.Rda")
FeaturePlot(Lung, features = c("CD68","FCGR3B","CD83",
                               "CD1C","HLA-DPA1","HLA-DPB1",
                               "KRT18","TPPP3","EPCAM",
                               "CD3D","CD3E",
                               "KLRD1",
                               "TPSB2",
                               "CD4","CD8A",
                               "MS4A1","CD79A","CD79B",
                               "IGHG4",                     
            label = T, order = T))
DotPlot(Lung, features=c("CD68","FCGR3B","CD83",
                         "CD1C","HLA-DPA1","HLA-DPB1",
                         "KRT18","TPPP3","EPCAM",
                         "CD3D","CD3E",
                         "KLRD1",
                         "TPSB2",
                         "CD4","CD8A",
                         "MS4A1","CD79A","CD79B",
                         "IGHG4",
                         
                         "CD4","CD8A","CD14","FCN1",
                         "MS4A1","CD79A","CD79B",
                         "IGHG4","JCHAIN"))
VlnPlot(Lung, features=c("CD68","FCGR3B","CD83",
                         "CD1C","HLA-DPA1","HLA-DPB1",
                         "KRT18","TPPP3","EPCAM",
                         "CD3D","CD3E",
                         "KLRD1",
                         "TPSB2",
                         "CD4","CD8A",
                         "MS4A1","CD79A","CD79B",
                         "IGHG4"),
        pt.size=0)
save(Lung,file="Lung2.Rda")
#自定义大群注释
new.cluster.ids <- c("Neutrophils","Macrophages","Neutrophils","Macrophages", #0-3          
                     "Neutrophils","Macrophages","T/NK","Macrophages", #4-7
                     "Macrophages","T/NK","Macrophages","Neutrophils","Macrophages", #8-12
                     "Neutrophils","Neutrophils","Macrophages","DCs",  #13-16
                     "Macrophages","Ciliated cells","Macrophages","Macrophages",#17-20
                     "Macrophages","Neutrophils","Macrophages","Macrophages","B cells",#21-25
                     "Secretory cells") 
names(new.cluster.ids) <- levels(Lung)
Lung <- RenameIdents(Lung, new.cluster.ids)
p2 <- DimPlot(Lung, reduction = "tsne", label = TRUE, pt.size = 0.5)
p1 <- DimPlot(Lung, reduction = "umap", label = TRUE, pt.size = 0.5)
p1 + p2
Lung[["celltype"]] <- Idents(Lung)
save(Lung, file = "Lung3.Rda")
