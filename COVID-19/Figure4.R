library(ggplot2)
library(Seurat)
library(monocle)
library(cowplot)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")
load("ISG-mono.Rda")
a1 <- plot_cell_trajectory(T_cds, color_by = "Pseudotime", cell_size=3,
                           cell_link_size=2)+facet_wrap("~group",
                                                                      nrow=1)
a1 <- a1 + theme(text=element_text(size=36),legend.text=element_text(size=20))
tiff("a1.tiff",width=40,height=30,units="cm",compression="lzw",res=300)
dev.new()
a1
dev.off()

a2 <- plot_cell_trajectory(T_cds, color_by = "celltype",cell_size=3,
                           cell_link_size=2)+facet_wrap("~group",
                                                                    nrow=1)
a2 <- a2 + theme(text=element_text(size=36),legend.text=element_text(size=20))
tiff("a2.tiff",width=40,height=30,units="cm",compression="lzw",res=300)
dev.new()
p2
dev.off()

##分支
BEAM_res <- BEAM(T_cds,branch_point=1,cores=8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name","pval","qval")]
BEAM_genes <- top_n(BEAM_res,n=50,desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_genes_branched_heatmap(T_cds[BEAM_genes,],
                                           branch_point=1,
                                           num_clusters=2,
                                           cores=8,
                                           use_gene_short_name=T,
                                           show_rownames=T,
                                  return_heatmap=T)
tiff("a3.tiff",width=40,height=30,units="cm",compression="lzw",res=300)
dev.new()
plot_heatmap(T_cds[BEAM_genes,],
                            branch_point=1,
                            num_clusters=2,
                            cores=8,
                            use_gene_short_name=T,
                            show_rownames=T,
                            return_heatmap=F)

dev.off()


##innate-like
load("T15.Rda")
a4 <- DimPlot(T15,label=T,label.size=10,pt.size=5)
a4 <- a4 + theme(text=element_text(size=36),legend.text=element_text(size=28))
tiff("a4.tiff",width=40,height=30,units="cm",compression="lzw",res=300)
dev.new()
a4
dev.off()

load("T15-mono.Rda")
a5 <- plot_cell_trajectory(T_cds, color_by = "celltype",cell_size=3,
                           cell_link_size=2)+facet_wrap("~group",nrow=1)
a5 <- a5 + theme(text=element_text(size=36),legend.text=element_text(size=20))
tiff(file='a5.tiff',compression="lzw",width=40,height=30,units="cm",res=300)
dev.new()
p8
dev.off()

#####
genes <- c("IFITM3","IFITM1")

a6 <- plot_genes_in_pseudotime(T_cds[genes,],color_by="group",cell_size=3)
a6 <- a6 + theme(text=element_text(size=36),legend.text=element_text(size=20))
tiff("a6.tiff",width=40,height=30,units="cm",compression="lzw",res=300)
dev.new()
a6
dev.off()
