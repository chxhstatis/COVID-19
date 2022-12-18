library(ggplot2)
library(Seurat)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/vaccine")
load("Bloodvac3.Rda")

#设置颜色
col<-c(brewer.pal(11,"Paired")[1:11])

f1 <- DimPlot(Blood,label=T,label.size=10,cols=col)
f1 <- f1 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f1.pdf",width=16,height=12)

##成分
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

##T_cell
load("T_cell2.Rda")
col2 <- col[3:9]

f3 <- DimPlot(T_cell,label=T,label.size=10,cols=col2)
f3 <- f3 + theme(text=element_text(size=36),legend.text=element_text(size=28))
ggsave("f3.pdf",width=16,height=12)

##成分
prop.table(table(Idents(T_cell)))
Cellratio <- prop.table(table(Idents(T_cell), T_cell$group), margin = 2)
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

