library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
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


#for (num in c(0.1,0.2,0.3,0.4)){
for (num in c(2)){
    pbmcnewer <- FindClusters(pbmcnewer, resolution = num)

    write.table(pbmcnewer@meta.data,paste('whole.stat.',num,'.for.raw.txt',sep=''),sep='\t',quote=F)

    #myresult<-Idents(pbmcnewer)
    myresult<-table(pbmcnewer$orig.ident,pbmcnewer@meta.data$seurat_clusters)
    write.table(myresult,file=paste('mergedata_clusters.r.',num,'.rmbatch.xls',sep=''),sep='\t',quote=F,row.names=T,col.names=T)

    p1 <- DimPlot(pbmcnewer, reduction = "umap", label = TRUE, repel = TRUE,raster=TRUE)
    oup<-paste('mergedata-UMAP-cluster-r',num,'.cluster.rmbatch.pdf',sep='')
    ggsave(oup, plot=p1, width = 8, height =8,limitsize=F)

	#DimPlot(pbmcnewer, reduction = "umap", label = TRUE, repel = TRUE)
    g1<-FeaturePlot(object = pbmcnewer, features = 'percent.mt', 
            cols = c("grey", "blue"), reduction = "umap",
            pt.size=1,min.cutoff=0,max.cutoff=5,ncol=1,raster=TRUE)
    ggsave('MT.feature.plot.pdf', plot=g1, width = 8, height = 8,limitsize=F)

    mymarkers_to_plot<-c('ENSMUSG00000053093','ENSMUSG00000057003','ENSMUSG00000033196','ENSMUSG00000056328','ENSMUSG00000028736','ENSMUSG00000020717','ENSMUSG00000019929')
    g1 <- FeaturePlot(object = pbmcnewer, features = mymarkers_to_plot, cols = c("grey", "blue"), reduction = "umap",pt.size=1,min.cutoff=0,ncol=4,max.cutoff=2,raster=TRUE)
    ggsave('All.gene.feature.plot.pdf', plot=g1, width = 32, height = 16,limitsize=F)
}


