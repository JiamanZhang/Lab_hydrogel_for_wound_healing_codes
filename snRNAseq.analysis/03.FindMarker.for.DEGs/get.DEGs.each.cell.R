library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
require("RColorBrewer")
library(dplyr)
library(pheatmap)

args<-commandArgs(TRUE)
cell<-args[1]
time1<-args[2]
time2<-args[3]


setwd('/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/02.snRNAseq.downstram.analysis/Rename.analysis/08.DEG.identified/Between.times.for.each.celltype.using.Seurat.for.G5/')
pbmcnewer<-readRDS('merge_processed.anchor.for.rmbatch.rename.rds')


pbmcnewer<-pbmcnewer[,pbmcnewer@meta.data$integrated_merge_cluster %in% c(cell)]

pbmcnewer<-pbmcnewer[,pbmcnewer@meta.data$treatment %in% c('G5')]

pbmcnewer<-pbmcnewer[,pbmcnewer@meta.data$time %in% c(time1,time2)]

dim(pbmcnewer@meta.data)


counts <- GetAssayData(object = pbmcnewer,assay='RNA', slot = "counts")
filtered_pbmc <- CreateSeuratObject(counts, meta.data = pbmcnewer@meta.data)

filtered_pbmc <- NormalizeData(filtered_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)



DEGdat<-FindMarkers(filtered_pbmc,assay='RNA',slot = 'counts',ident.1=time1,ident.2=time2,group.by='time',logfc.threshold =0,min.pct = 0.1)


oup<-paste('./min.pct.0.1/',cell,'.',time1,'.',time2,'.DEGs.result.G5.txt',sep='')
write.table(DEGdat,oup,quote=F,sep='\t')

