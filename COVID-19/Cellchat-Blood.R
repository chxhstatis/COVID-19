library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(svglite)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")
load("Blood3.Rda")
cell <- subset(Blood,group=='control')
##control
##提取表达矩阵和细胞分类信息
cell@meta.data[["celltype"]]<- as.factor(as.character(cell@meta.data[["celltype"]]))
#创建cellchat对象
cellchat <- createCellChat(object = cell@assays[["RNA"]]@data, meta = cell@meta.data, group.by = "celltype") 
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling") 
#预处理表达数据以进行细胞间相互作用分析

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat<- projectData(cellchat,PPI.human )
#推断细胞间相互作用网络与分析

#根据表达值推测细胞互作的概率（cellphonedb哦嗯平均值代表互作强度）
cellchat <- computeCommunProb(cellchat,raw.use = F,population.size = T)
#如果不想用上一步PPI的结果，raw.use = T
cellchat<- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat<- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

par(mfrow=c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight=groupSize,weight.scale=T,
                 label.edge=F,title.name="Number of interactions",vertex.label.cex=2)
netVisual_circle(cellchat@net$weight, vertex.weight=groupSize,weight.scale=T,
                 label.edge=F,title.name="Interaction Weights",vertex.label.cex=2)
cellchat@netP$pathways
levels(cellchat@idents)

cellchat <- netAnalysis_computeCentrality(cellchat,slot.name="netP")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

ht1 <- netAnalysis_signalingRole_heatmap(cellchat,pattern="outgoing",height=14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,pattern="incoming",height=14)
ht1+ht2

nPatterns = 5
myidentifyCommunicationPatterns <- edit(identifyCommunicationPatterns)
environment(myidentifyCommunicationPatterns) <- environment(myidentifyCommunicationPatterns)
cellchat <- myidentifyCommunicationPatterns(cellchat, pattern="outgoing",
                                            k = nPatterns,height=12)
netAnalysis_river(cellchat,pattern="outgoing")
netAnalysis_dot(cellchat,pattern="outgoing")

cellchat <- myidentifyCommunicationPatterns(cellchat, pattern="incoming",
                                            k = nPatterns,height=12)
netAnalysis_river(cellchat,pattern="incoming")
netAnalysis_dot(cellchat,pattern="incoming")

save(cellchat,file='cellchat_c.Rda')

####patient
cell <- subset(Blood,group!='control')

##提取表达矩阵和细胞分类信息
cell@meta.data[["celltype"]]<- as.factor(as.character(cell@meta.data[["celltype"]]))
#创建cellchat对象
cellchat <- createCellChat(object = cell@assays[["RNA"]]@data, meta = cell@meta.data, group.by = "celltype") 
cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling") 
#预处理表达数据以进行细胞间相互作用分析

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat<- projectData(cellchat,PPI.human )
#推断细胞间相互作用网络与分析

#根据表达值推测细胞互作的概率（cellphonedb哦嗯平均值代表互作强度）
cellchat <- computeCommunProb(cellchat,raw.use = F,population.size = T)
#如果不想用上一步PPI的结果，raw.use = T
cellchat<- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat<- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

par(mfrow=c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight=groupSize,weight.scale=T,
                 label.edge=F,title.name="Number of interactions",vertex.label.cex=2)
netVisual_circle(cellchat@net$weight, vertex.weight=groupSize,weight.scale=T,
                 label.edge=F,title.name="Interaction Weights",vertex.label.cex=2)
cellchat@netP$pathways
levels(cellchat@idents)

cellchat <- netAnalysis_computeCentrality(cellchat,slot.name="netP")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1

ht1 <- netAnalysis_signalingRole_heatmap(cellchat,pattern="outgoing",height=14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,pattern="incoming",height=14)
ht1+ht2

nPatterns = 5
myidentifyCommunicationPatterns <- edit(identifyCommunicationPatterns)
environment(myidentifyCommunicationPatterns) <- environment(myidentifyCommunicationPatterns)
cellchat <- myidentifyCommunicationPatterns(cellchat, pattern="outgoing",
                                            k = nPatterns,height=12)
netAnalysis_river(cellchat,pattern="outgoing")
netAnalysis_dot(cellchat,pattern="outgoing")

cellchat <- myidentifyCommunicationPatterns(cellchat, pattern="incoming",
                                            k = nPatterns,height=12)
netAnalysis_river(cellchat,pattern="incoming")
netAnalysis_dot(cellchat,pattern="incoming")

save(cellchat,file='cellchat_p.Rda')