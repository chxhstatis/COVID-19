library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)
library(reshape2)

setwd("/run/media/xuhao/b218f49e-a30d-4700-a6bf-e8186e6c9310/nf/COVID-19_single-cell/Blood/patient/predata")
load('T_cell2.Rda')
ISG <- subset(T_cell,celltype=='ISG CD4T'|celltype=='ISG CD8T'|celltype=='innate-like Tcells')
data <- ISG@assays$RNA@counts
pd <- new('AnnotatedDataFrame', data = ISG@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
T_cds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
#--------------------------------------------------------------
#评估离散度
T_cds <- estimateSizeFactors(T_cds)
T_cds <- estimateDispersions(T_cds)
#过滤
T_cds <- detectGenes(T_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(T_cds),
                                    num_cells_expressed >= 10))
length(expressed_genes)
#分布
L <- log(exprs(T_cds[expressed_genes,]))
M <- Matrix::t(scale(Matrix::t(L)))
melted_dens_df <- melt(M)
qplot(value,geom="density",data=melted_dens_df)+
  stat_function(fun=dnorm,size=0.5,color="red")+
  xlab("Standardized log(FPKM)")+
  ylab("Density")
disp_table <- dispersionTable(T_cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
T_cds <- setOrderingFilter(T_cds, disp.genes)
plot_ordering_genes(T_cds)
#降维
T_cds <- reduceDimension(T_cds, max_components = 2,
                         method = "DDRTree")
#排列
T_cds <- orderCells(T_cds)
#按“Pseudotime”分组；
plot_cell_trajectory(T_cds, color_by = "IFITM3")
plot_cell_trajectory(T_cds, color_by = "Pseudotime")+facet_wrap("~group",
                                                                nrow=1)
#按“State”分组；
plot_cell_trajectory(T_cds, color_by = "State")
#按seurat分群结果分组
plot_cell_trajectory(T_cds, color_by = "seurat_clusters")
#按细胞类型分组
plot_cell_trajectory(T_cds, color_by = "celltype")+facet_wrap("~group",
                                                              nrow=1)

plot_cell_trajectory(T_cds, color_by="State")+facet_wrap("~group",
                                                            nrow=1)
plot_cell_trajectory(T_cds, color_by="State")+facet_wrap("~celltype",
                                                         nrow=1)
plot_cell_trajectory(T_cds, color_by="group")+facet_wrap("~seurat_clusters",
                                                         nrow=1)
plot_cell_trajectory(T_cds, color_by="seurat_clusters")+facet_wrap("~group",
                                                         nrow=1)
de <- differentialGeneTest(T_cds, fullModelFormulaStr=" ~sm.ns(Pseudotime)",cores=8)
de_genes <- top_n(de,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_pseudotime_heatmap(T_cds[de_genes,],
                        num_clusters=3,
                        cores=8,
                        show_rownames=T)

save(T_cds, file = "ISG-mono.Rda")


#####T15
load('T15.Rda')
data <- T15@assays$RNA@counts
pd <- new('AnnotatedDataFrame', data = T15@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
T_cds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

#评估离散度
T_cds <- estimateSizeFactors(T_cds)
T_cds <- estimateDispersions(T_cds)
#过滤
T_cds <- detectGenes(T_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(T_cds),
                                    num_cells_expressed >= 10))
length(expressed_genes)
#分布
L <- log(exprs(T_cds[expressed_genes,]))
M <- Matrix::t(scale(Matrix::t(L)))
melted_dens_df <- melt(M)
qplot(value,geom="density",data=melted_dens_df)+
  stat_function(fun=dnorm,size=0.5,color="red")+
  xlab("Standardized log(FPKM)")+
  ylab("Density")
disp_table <- dispersionTable(T_cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
T_cds <- setOrderingFilter(T_cds, disp.genes)
plot_ordering_genes(T_cds)
#降维
T_cds <- reduceDimension(T_cds, max_components = 2,
                         method = "DDRTree")
#排列
T_cds <- orderCells(T_cds)
#按“Pseudotime”分组；
plot_cell_trajectory(T_cds, color_by = "IFITM3")
plot_cell_trajectory(T_cds, color_by = "Pseudotime")+facet_wrap("~group",
                                                                nrow=1)
#按“State”分组；
plot_cell_trajectory(T_cds, color_by = "State")
#按seurat分群结果分组
plot_cell_trajectory(T_cds, color_by = "seurat_clusters")
#按细胞类型分组
plot_cell_trajectory(T_cds, color_by = "celltype")+facet_wrap("~group",
                                                              nrow=1)

plot_cell_trajectory(T_cds, color_by="State")+facet_wrap("~group",
                                                         nrow=1)
plot_cell_trajectory(T_cds, color_by="State")+facet_wrap("~celltype",
                                                         nrow=1)
plot_cell_trajectory(T_cds, color_by="group")+facet_wrap("~seurat_clusters",
                                                         nrow=1)
plot_cell_trajectory(T_cds, color_by="seurat_clusters")+facet_wrap("~group",
                                                                   nrow=1)
de <- differentialGeneTest(T_cds, fullModelFormulaStr=" ~sm.ns(Pseudotime)",cores=8)
de_genes <- top_n(de,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
plot_pseudotime_heatmap(T_cds[de_genes,],
                        num_clusters=3,
                        cores=8,
                        show_rownames=T)

save(T_cds, file = "T15-mono.Rda")