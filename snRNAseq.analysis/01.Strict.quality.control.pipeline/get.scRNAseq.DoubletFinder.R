library(DoubletFinder)
# library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)

args<-commandArgs(TRUE)
sample<-args[1]
input<-args[2]

setwd(paste(input,'/DoubletFinder.result/',sample,'/',sep=''))

tocinp<-paste(input,'/soupX_matrix/',sample,'/matrix/',sep='')

toc <- Read10X(tocinp,gene.column=2)

all <- toc
all <- CreateSeuratObject(all)
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)
all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:30)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:30)

keloid<-all

sweep.res.list_keloid <- paramSweep_v3(keloid, PCs = 1:30, sct = FALSE)
sweep.stats_keloid <- summarizeSweep(sweep.res.list_keloid, GT = FALSE)
bcmvn_keloid <- find.pK(sweep.stats_keloid) 

mpK<-as.numeric(as.vector(bcmvn_keloid$pK[which.max(bcmvn_keloid$BCmetric)]))

DoubletRate = ncol(keloid)*8*1e-6

annotations <- keloid@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  

nExp_poi <- round(DoubletRate*length(keloid$seurat_clusters))

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

keloid <- doubletFinder_v3(keloid, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
keloid <- doubletFinder_v3(keloid, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)

keloid@meta.data[,"DF.nExp_poi"]<-keloid@meta.data[,paste('DF.classifications_0.25_',mpK,'_',nExp_poi,sep='')]
keloid@meta.data[,"DF.nExp_poi.adj"]<-keloid@meta.data[,paste('DF.classifications_0.25_',mpK,'_',nExp_poi.adj,sep='')]

# keloid@meta.data[,"DF_hi.lo"] <- keloid@meta.data$DF.classifications_0.25_0.02_184
keloid@meta.data[,"DF_hi.lo"]<-keloid@meta.data$DF.nExp_poi

keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet" & keloid@meta.data$DF.nExp_poi.adj == "Singlet")] <- "Doublet-Low Confidience"
keloid@meta.data$DF_hi.lo[which(keloid@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
table(keloid@meta.data$DF_hi.lo)

#oud<-paste(input,'/',sample,'/',sep='')

p1<-DimPlot(keloid, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))
ggsave(paste(sample,'.Doublet.pdf',sep=''), plot=p1, width = 8, height =6)

p1 <- DimPlot(keloid, reduction = "umap", label = TRUE, repel = TRUE)
ggsave(paste(sample,'.UMAP.pdf',sep=''), plot=p1, width = 7, height =6)

write.table(keloid@meta.data,paste(sample,'.Doublet.info.txt',sep=''),sep='\t',quote=F)


