library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(cowplot)
#library(limma)
library(harmony)
#library(DESeq2)
#library(scater)
#library(SingleR)
#library(scRNAseq)
library(dplyr)
#library(pheatmap)
#library(SeuratData)

#inp<-'/Lustre03/data/zhangJM/01.Projects/16.Pig.dif.Muscle.snRNA.for.JIN/test/Exp.matrix/whole.new.samples.txt'
#sdat<-read.table(inp,head=F,sep='\t',stringsAsFactors=F)
setwd('/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/01.snRNAseq.integrate.pipeline/')

# sample_list<-c('BY2_GM1','BY4_GM1','LM1_GM1','LM2_GM1','BY2_PM1','BY4_PM1','LM1_PM1','LM2_PM1')
inp<-'samples.txt'
dat<-read.table(inp,head=F,sep='\t')
sample_list<-dat$V1

dir<-c()
for (sample in sample_list){
    ind<-'/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/01.snRNAseq.integrate.pipeline/soupX_matrix/'
    dirname<-paste(ind,sample,'/matrix/',sep='')
    dir<-c(dir,dirname)
}


names(dir) = gsub('_','-',sample_list)
scRNAlist <- list()
for (i in 1:length(sample_list)){
	counts <- Read10X(data.dir = dir[i])
	# scRNAlist[[i]] <- CreateSeuratObject(counts, project =sample_list[i], min.cells = 3, min.features = 10)
    scRNAlist[[i]] <- CreateSeuratObject(counts, project =sample_list[i])
    head(scRNAlist[[i]]@meta.data, 2)

    inp<-paste('/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/01.snRNAseq.integrate.pipeline/DoubletFinder.result/',sample_list[i],'/',sample_list[i],'.Doublet.info.txt',sep='')
    DFdat<-read.table(inp,head=T,sep='\t',stringsAsFactors = F)
    DFcell_list<-paste(gsub('_','-',sample_list[i]),'_',rownames(DFdat[DFdat$DF_hi.lo=='Singlet',]),sep='')

    # scRNAlist[[i]]

    scRNAlist[[i]]<-subset(scRNAlist[[i]],cells=c(DFcell_list))

    # scRNAlist[[i]]


    Mgid_list<-c('ENSMUSG00000064341','ENSMUSG00000064345','ENSMUSG00000064351','ENSMUSG00000064354','ENSMUSG00000064356','ENSMUSG00000064357','ENSMUSG00000064358','ENSMUSG00000064360','ENSMUSG00000065947','ENSMUSG00000064363','ENSMUSG00000064367','ENSMUSG00000064368','ENSMUSG00000064370')

    num<-0
    for (Mgid in Mgid_list){
        if (Mgid %in% rownames(scRNAlist[[i]]@assays$RNA)) {
            num=num+1
            scRNAlist[[i]][[Mgid]] <- PercentageFeatureSet(scRNAlist[[i]], assay = "RNA", feature = Mgid)
        }
    }

    scRNAlist[[i]][["percent.mt"]] <- apply(scRNAlist[[i]]@meta.data[,4:(4+num-1)], 1, sum)
    scRNAlist[[i]]@meta.data <- scRNAlist[[i]]@meta.data[,-c(4:(4+num-1))]

    print(sample_list[i])
    print(num)
    
    table(scRNAlist[[i]]@meta.data$orig.ident)

}

scRNA2 <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]],scRNAlist[[3]],scRNAlist[[4]],scRNAlist[[5]],scRNAlist[[6]],scRNAlist[[7]],scRNAlist[[8]],scRNAlist[[9]],scRNAlist[[10]],scRNAlist[[11]],scRNAlist[[12]]),add.cell.ids = sample_list,project='pbmc3k')


table(scRNA2@meta.data$orig.ident)
# write.table(scRNA2@meta.data$orig.ident,'raw.cell.num.txt',quote=F,sep='\t')

head(scRNA2@meta.data, 2)


rb.genes <- rownames(scRNA2)[grep("^RP[SL]",rownames(scRNA2))]
scRNA2[["percent.rb"]] <- PercentageFeatureSet(scRNA2, pattern="^RP[SL]")


#scRNA2<-subset(scRNA2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 5)

pbmc<-scRNA2
table(pbmc@meta.data$orig.ident)
head(pbmc@meta.data, 2)

pbmcnew<-pbmc

inp<-'Mm39.110.PCG.no.MT.genes.txt'
dat<-read.table(inp,head=F,sep='\t')

MT_list<-dat$V4

pbmcnew
counts <- GetAssayData(object = pbmcnew, slot = "counts")
filter_counts<-counts[rownames(counts)%in%MT_list,]
filtered_pbmc <- CreateSeuratObject(filter_counts, meta.data = pbmcnew@meta.data)

filtered_pbmc
print('1')

bm280k<-filtered_pbmc
bm280k
bm280k <- SCTransform(bm280k,assay = "RNA")
bm280k <- RunPCA(bm280k,assay='SCT', verbose = FALSE)
bm280k.integrated<-RunHarmony(bm280k,group.by.vars = 'orig.ident',assay.use='SCT',max.iter.harmony=30)
bm280k.integrated <- RunTSNE(bm280k.integrated,reduction='harmony', dims = 1:30)
bm280k.integrated <- RunUMAP(bm280k.integrated, reduction = "harmony", dims = 1:30)


setwd('/Lustre03/data/zhangJM/01.Projects/20.Mouse.SMT.snRNAseq.for.Haohuan.Li/01.snRNAseq.integrate.pipeline/raw.UMAP.for.MT.filter/')

plot1 <- DimPlot(bm280k.integrated, reduction = "tsne", group.by="orig.ident",raster=TRUE)
ggsave('mergedata.TSNE.cluster.all.raw.pdf', plot=plot1, width = 8, height =8,limitsize=F)

plot1 <- DimPlot(bm280k.integrated, reduction = "tsne", split.by="orig.ident",group.by="orig.ident",ncol=6,raster=TRUE)
ggsave('mergedata.TSNE.cluster.all.split.samples.raw.pdf', plot=plot1, width = 48, height =16,limitsize=F)


plot1 <- DimPlot(bm280k.integrated, reduction = "umap", group.by="orig.ident",raster=TRUE)
ggsave('mergedata.UMAP.cluster.all.raw.pdf', plot=plot1, width = 8, height =8,limitsize=F)

plot1 <- DimPlot(bm280k.integrated, reduction = "umap", split.by="orig.ident",group.by="orig.ident",ncol=6,raster=TRUE)
ggsave('mergedata.UMAP.cluster.all.split.samples.raw.pdf', plot=plot1, width = 48, height =16,limitsize=F)


saveRDS(bm280k.integrated, file = "merge_processed.anchor.for.rmbatch.raw.cells.rds")


