library(ggplot2)
library(Seurat)
library(DAseq)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")
load("T_cell2.Rda")

##设置颜色
col<-c(brewer.pal(11,"Paired")[1:15])
col[12] <- '#9933FA'
col[13] <- '#FA8072'
col[14] <- '#FF6103'
col[15] <- '#03A89E'


p1 <- DimPlot(T_cell,label=T,label.size=10,cols=col)
p1 <- f1 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("p1.pdf",width=16,height=12)

##成分
prop.table(table(Idents(T_cell)))
Cellratio <- prop.table(table(Idents(T_cell), T_cell$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","group","Freq")
colourCount = length(unique(Cellratio$celltype))
p2 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Freq, y= group, fill = celltype),stat = "identity",width = 0.5,size = 0.1,colour = '#000000')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Group')+
  coord_flip()
f2 <- f2 + scale_fill_manual(values=col)
p2 <- p2 + theme(text=element_text(size=36))
ggsave("p2.pdf",width=16,height=12)


###control/patient
T_cell$group[T_cell$group=='mild-moderate'] <- 'patient'
T_cell$group[T_cell$group=='severe'] <- 'patient'
table(T_cell@meta.data$group)
#设置一些基本的 Python 参数
python2use <- "/usr/bin/python3.6"
GPU <- 5

DimPlot(T_cell, reduction = "umap",group.by='celltype',label=T,cols=col)
#获取标签信息
label2<-cbind(colnames(T_cell), T_cell$group)
label2<-as.data.frame(label2)
colnames(label2)<-c("cells","group")
labels_res <- label2[label2$group == "control", "cells"]
labels_nonres <- label2[label2$group ==  "patient", "cells"]

#获取 DA 细胞
dat1<-T_cell@reductions$harmony@cell.embeddings
dat1<-as.data.frame(dat1)
dat.umap<-T_cell@reductions$umap@cell.embeddings

da_cells <- getDAcells(
  X = dat1,
  cell.labels = label2$cells,
  labels.1 = labels_res,
  labels.2 = labels_nonres,
  k.vector = seq(50, 500, 50),
  plot.embedding = dat.umap
)
str(da_cells[1:2])
da_cells$pred.plot
da_cells$rand.plot

da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(-0.93,0.94),
  plot.embedding = dat.umap
)
p3 <- da_cells$da.cells.plot
p3 <- p3 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("p3.pdf",width=16,height=12)

#获取 DA 区域
da_regions <- getDAregion(
  X = dat1,
  da.cells = da_cells,
  cell.labels =label2$cells,
  labels.1 = labels_res,
  labels.2 = labels_nonres,
  resolution = 0.01,
  plot.embedding =dat.umap,
)
str(da_regions[1:2])
a <- da_regions$DA.stat
da_regions$da.region.plot

T_cell@meta.data[["da.cluster"]]<-da_regions[["da.region.label"]]
T_cell@meta.data[["da.scores"]]<-dascores
table(T_cell@meta.data[["da.cluster"]])
DimPlot(T_cell,group.by="da.cluster",label=T)


table(da_regions[["da.region.label"]])
T_cell@meta.data[["da.cluster"]]<-da_regions[["da.region.label"]]
table(T_cell@meta.data[["da.cluster"]])

col2 <- c(brewer.pal(5,"Paired")[1:5])
col2[1] <- '#C0C0C0'
col2[4] <- '#B0171F'
col2[5] <- '#FF7F50'

p4 <- DimPlot(T_cell, group.by="da.cluster",label=T,cols=col2,label.size=20)
p4 <- p4 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("p4.pdf",width=16,height=12)

#marker
a <- T_cell@assays$RNA@data
STG_marker <- STGmarkerFinder(
  X = a,
  da.regions=da_regions,
  lambda=1.5, n.runs=5, return.model=T,
  python.use=python2use, GPU=GPU
)
save(STG_marker,file="marker.Rda")
p5 <- DotPlot(T_cell,group.by='da.cluster',features=c("LTB","CD69","IRF7","NKG7","CCL5","CX3CR1","FCER1G","FCGR3B","S100A8",
                                                "IFITM1","IFITM2","IFITM3","SELL","CCR7","TCF7"))
p5 <- p5 + ylab("DA subpopulations")+xlab("")
p5 <- p5 +theme(text=element_text(size=36),legend.text=element_text(size=12))
p5 <- p5 +theme(axis.text.x = element_text(angle = 0, vjust = 1, 
                                           size = 28, hjust = 0.5),
                axis.text.y = element_text(angle = 0, vjust = 1, 
                                           size = 34, hjust = 0.5))
ggsave("p5.pdf",width=30,height=12)
