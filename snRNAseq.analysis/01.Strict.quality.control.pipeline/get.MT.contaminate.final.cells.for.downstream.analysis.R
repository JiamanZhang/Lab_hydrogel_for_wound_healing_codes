library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
library(DropletUtils)
# library(limma)
#library(DESeq2)
#library(scater)
#library(SingleR)
#library(scRNAseq)
library(dplyr)
#library(pheatmap)
#library(SeuratData)

setwd('/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/01.snRNAseq.integrate.pipeline/raw.UMAP.for.MT.filter/')

pbmcnewer<-readRDS('merge_processed.anchor.for.rmbatch.raw.cells.rds')
pbmcnewer <- FindNeighbors(pbmcnewer, reduction = "harmony", dims = 1:30)

num<-2
pbmcnewer <- FindClusters(pbmcnewer, resolution = num)

pbmc_small<-pbmcnewer

#pbmc_small<-subset(pbmc_small,idents=c('6'),invert=T)

pbmc_small <- RenameIdents(pbmc_small, '0'='others','1'='others','2'='others','3'='others','4'='others','5'='others','6'='others','7'='others','8'='others','9'='others','10'='others','11'='others','12'='others','13'='others','14'='others','15'='others','16'='others','17'='others','18'='others','19'='others','20'='others','21'='others','22'='others','23'='others','24'='others','25'='others','26'='others','27'='others','28'='others','29'='others','30'='MT-30','31'='others','32'='others','33'='others','34'='others','35'='others','36'='others','37'='others','38'='others','39'='others','40'='others','41'='MT-41','42'='others','43'='others','44'='others','45'='others')
levels(pbmc_small)

order_list<-c('MT-30','MT-41','others')
pbmc_small$integrated_merge_cluster<-factor(x =Idents(pbmc_small), levels = order_list)

p1 <- DimPlot(pbmc_small, reduction = "umap",group.by='integrated_merge_cluster', label = TRUE, repel = TRUE,raster=TRUE)
oup<-paste('mergedata-UMAP-cluster-r',num,'.cluster.rmbatch.rename.pdf',sep='')
ggsave(oup, plot=p1, width = 10, height =8,limitsize=F)


pbmcnewer<-subset(pbmcnewer,idents=c('30','41'),invert=T)

g1<-FeaturePlot(object = pbmcnewer, features = 'percent.mt', 
            cols = c("grey", "blue"), reduction = "umap",
            pt.size=1,min.cutoff=0,max.cutoff=5,ncol=1,raster=TRUE)
ggsave(paste('mergedata-UMAP-cluster-r.',num,'.cluster.rmbatch.percent.mt.filter.MT.cluster.pdf',sep=''), plot=g1, width = 8, height = 8,limitsize=F)

p1 <- DimPlot(pbmcnewer, reduction = "umap", label = TRUE, repel = TRUE,raster=TRUE)
oup<-paste('mergedata-UMAP-cluster-r',num,'.cluster.rmbatch.filter.MT.cluster.pdf',sep='')
ggsave(oup, plot=p1, width = 8, height =8,limitsize=F)


write.table(pbmcnewer@meta.data,'merge_processed.anchor.for.rmbatch.raw.cells.filter.MT.meta.data.txt',quote=F,sep='\t')

saveRDS(pbmcnewer, file = "merge_processed.anchor.for.rmbatch.raw.cells.filter.MT.cluster.rds")

oup<-paste('/Lustre02/zhangJM/01.Project/16.Mouse.SMT.snRNAseq.for.Haohuan.Li/01.Mouse.muscle.snRNAseq.analysis/01.snRNAseq.integrate.pipeline/raw.UMAP.for.MT.filter/filt.MTcluster.matrix/',sep='')
write10xCounts(oup, pbmcnewer@assays$RNA@counts,version="3")

