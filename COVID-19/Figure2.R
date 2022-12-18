library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(CellChat)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")
load("Blood3.Rda")

#设置颜色
col<-c(brewer.pal(11,"Paired")[1:13])
col[12] <- '#9933FA'
col[13] <- '#F5F5F5'

f1 <- DimPlot(Blood,label=T,label.size=14,cols=col)
f1 <- f1 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f1.pdf",width=16,height=12)

#成分
prop.table(table(Idents(Blood)))
Cellratio <- prop.table(table(Idents(Blood), Blood$group), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","group","Freq")
colourCount = length(unique(Cellratio$celltype))
f2 <- ggplot(Cellratio) + 
  geom_bar(aes(x =Freq, y= group, fill = celltype),stat = "identity",width = 0.5,size = 0.1,colour = '#000000')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Group')+
  coord_flip()
f2 <- f2 + scale_fill_manual(values=col)
f2 <- f2 + theme(text=element_text(size=36))
ggsave("f2.pdf",width=16,height=12)

##--------------------------------------------------------------
#通讯
load("cellchat_c.Rda")
cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name="netP")

f3 <- netAnalysis_signalingRole_scatter(cellchat,font.size=24, dot.size=(c(5,20)),label.size=20)
ggsave("f3.pdf",width=16,height=12)

f5 <- netVisual_aggregate(cellchat, signaling = "IL6",vertex.label.cex=2, layout = "circle",pt.title=30)
f7 <- netVisual_aggregate(cellchat, signaling = "TNF",vertex.label.cex=2, layout = "circle",pt.title=30)

###patient
load("cellchat_P.Rda")
cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name="netP")

f4 <- netAnalysis_signalingRole_scatter(cellchat,font.size=24, dot.size=(c(5,20)),label.size=20)
ggsave("f4.pdf",width=16,height=12)

f6 <- netVisual_aggregate(cellchat, signaling = "IL6",vertex.label.cex=2, layout = "circle",pt.title=30)
f8 <- netVisual_aggregate(cellchat, signaling = "TNF",vertex.label.cex=2, layout = "circle",pt.title=30)

