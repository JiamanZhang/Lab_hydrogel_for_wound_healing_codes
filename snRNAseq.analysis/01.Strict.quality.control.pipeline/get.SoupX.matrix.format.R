library(SoupX)
library(Seurat)
library(DropletUtils)

args<-commandArgs(TRUE)
sample<-args[1]
input<-args[2]

setwd(paste(input,'/soupX_matrix/',sample,'/',sep=''))

ind<-paste(input,'/matrix/',sample,'/',sep='')

tocinp<-paste(ind,'/filtered.matrix/',sep='')
todinp<-paste(ind,'/raw.matrix/',sep='')

toc <- Read10X(tocinp,gene.column=1)
#table(toc)
tod <- Read10X(todinp,gene.column=1)
tod <- tod[rownames(toc),]

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

matx <- all@meta.data
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)

rhonum<-unique(sc$metaData$rho)
sc = setContaminationFraction(sc, rhonum)
out = adjustCounts(sc)

oup<-paste(input,'/soupX_matrix/',sample,'/matrix/',sep='')
DropletUtils:::write10xCounts(oup, out,version="3")

rhonum_dat<-data.frame(sample=sample,rho=rhonum)
oup<-paste(input,'/soupX_matrix/',sample,'/',sample,'.SoupX.rho.txt',sep='')
write.table(rhonum_dat,oup,quote=F,sep='\t',row.names=F)



