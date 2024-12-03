library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
require("RColorBrewer")
library(dplyr)
library(pheatmap)
args<-commandArgs(TRUE)
treatment<-'G5'
orig.ident<-args[1]
celltype<-args[2]



setwd('/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/02.snRNAseq.downstram.analysis/Rename.analysis/07.pseudoBulk.for.each.cell.type/')

pbmcnewer<-readRDS('merge_processed.anchor.for.rmbatch.rename.rds')
pbmcnewer<-pbmcnewer[,pbmcnewer@meta.data$treatment %in% c('G5')]

newpbmcnewer<-pbmcnewer[,pbmcnewer@meta.data$orig.ident==orig.ident&pbmcnewer@meta.data$integrated_merge_cluster==celltype]

#oup<-paste('./result.data/',tissue,'.',breed,'.',celltype,'.snRNAseq.meta.data.cell.info.txt',sep='')
#write.table(newpbmcnewer@meta.data,oup,sep='\t',quote=F)

#oup<-paste('./result.data/',tissue,'.',breed,'.',celltype,'.whole.RNA.gene.counts.csv',sep='')
#write.csv(file=oup,as.data.frame(newpbmcnewer[["RNA"]]@counts),quote = F)

oup<-paste('./result.data/',orig.ident,'.',celltype,'.whole.RNA.gene.data.csv',sep='')
write.csv(file=oup,as.data.frame(newpbmcnewer[["RNA"]]@data),quote = F)

