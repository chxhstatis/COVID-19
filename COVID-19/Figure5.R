library(ggplot2)
library(Seurat)
library(DAseq)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Lung/patient/predata")
load("Lung3.Rda")

##设置颜色
col<-c(brewer.pal(7,"Paired")[1:7])

##----
f1 <- DimPlot(Lung,label=T,label.size=10,cols=col)
f1 <- f1 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f1.pdf",width=16,height=12)

##成分
prop.table(table(Idents(Lung)))
Cellratio <- prop.table(table(Idents(Lung), Lung$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","group","Freq")
colourCount = length(unique(Cellratio$celltype))
f2 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Freq, y= group, fill = celltype),stat = "identity",width = 0.7,size = 0.3,colour = '#222222')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Group')+
  coord_flip()
f2 <- f2 + theme(text=element_text(size=22))
f2 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Freq, y= group, fill = celltype),stat = "identity",width = 0.5,size = 0.1,colour = '#000000')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Group')+
  coord_flip()
f2 <- f2 + scale_fill_manual(values=col)
f2 <- f2 + theme(text=element_text(size=36))
ggsave("f2.pdf",width=16,height=12)

##T/NK
load("T_NKnew3.Rda")
col2 <- c(brewer.pal(9,"Paired")[2:10])

f3 <- DimPlot(T_NK,label=T,label.size=10,cols=col2)
f3 <- f3 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f3.pdf",width=16,height=12)

##成分
prop.table(table(Idents(T_NK)))
Cellratio <- prop.table(table(Idents(T_NK), T_NK$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","group","Freq")
colourCount = length(unique(Cellratio$celltype))
f4 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Freq, y= group, fill = celltype),stat = "identity",width = 0.5,size = 0.1,colour = '#000000')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Group')+
  coord_flip()
f4 <- f4 + scale_fill_manual(values=col2)
f4 <- f4 + theme(text=element_text(size=36))
ggsave("f4.pdf",width=16,height=12)

###control/patient
table(T_NK@meta.data$group)
#设置一些基本的 Python 参数
python2use <- "/usr/bin/python3.6"
GPU <- 5

DimPlot(T_NK, reduction = "umap",group.by='celltype',label=T)
#获取标签信息
label2<-cbind(colnames(T_NK), T_NK$group)
label2<-as.data.frame(label2)
colnames(label2)<-c("cells","group")
labels_res <- label2[label2$group == "control", "cells"]
labels_nonres <- label2[label2$group ==  "patient", "cells"]

#获取 DA 细胞
dat1<-T_NK@reductions$harmony@cell.embeddings
dat1<-as.data.frame(dat1)
dat.umap<-T_NK@reductions$umap@cell.embeddings

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
  X = da_cells, pred.thres = c(-0.8,0.86),
  plot.embedding = dat.umap
)
f5 <- da_cells$da.cells.plot
f5 <- f5+theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f5.pdf",width=16,height=12)
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
da_regions$DA.stat
da_regions$da.region.plot

col2[1] <- '#DCDCDC'
col2[2] <- '#FF3030'
col2[4] <- '#0000CD'
T_NK@meta.data[["da.cluster"]]<-da_regions[["da.region.label"]]
table(T_NK@meta.data[["da.cluster"]])
f6 <- DimPlot(T_NK, group.by="da.cluster",label=T,cols=col2,label.size=20)
f6 <- f6 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f6.pdf",width=16,height=12)

#marker
a <- T_cell@assays$RNA@data
STG_marker <- STGmarkerFinder(
  X = a,
  da.regions=da_regions,
  lambda=1.5, n.runs=5, return.model=T,
  python.use=python2use, GPU=GPU
)
save(STG_marker,file="T/NK_marker.Rda")

