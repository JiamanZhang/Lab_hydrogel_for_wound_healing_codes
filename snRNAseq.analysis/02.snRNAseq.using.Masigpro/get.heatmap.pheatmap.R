library(reshape2)
require("RColorBrewer")
library(gplots)
library(pheatmap)
library(reshape2)
library(plyr)
library(ape)
library(scales)
library(gtools)
library(RColorBrewer)
library(preprocessCore)
setwd('/Users/zhangjiaman/Desktop/students.works/zjm.data/hic项目/34.Mouse.LDM.snRNAseq.for.HaoHuan/New.snRNAseq.analysis/02.snRNAseq.downstram.analysis/Rename/07.pseudoBulk.for.each.cell.type/01.Masigpro/Immune/')

inp<-'Immune.TPM.cluster.sorted.txt'
dat<-read.table(inp,head=T,sep='\t',stringsAsFactors = F)
#dat<-dat[dat$gid%in%gdat$V2,]
order_list<-dat$cluster
rownames(dat)<-dat$gid
#rownames(dat)<-gname_list
dat<-dat[,-c(1,8)]

annotation_row = data.frame(order=order_list)
rownames(annotation_row) <- rownames(dat)

annotation_col = data.frame(time=c(rep('D3',time=2),rep('D9',time=2),rep('D13',time=2)))
rownames(annotation_col) <- colnames(dat)


annotation_colors<-list(order=c('1'='#009000','2'='#B362F1','3'='#FF9E04','4'='#3179F0'),
                        time=c(D3='#0095FF',D9='#FA6DDD',D13='#F1987A'
                        ))


break_list<-seq(-1,1,0.01)
pdf('test.heatmap.pdf',w=4.5,h=6)
pheatmap(dat,cluster_rows = F,color=colorRampPalette(c('#2181c5','white','#cc261c'))(length(break_list)),
         cluster_cols=F,show_rownames = F,show_colnames=T,breaks=break_list,
         #annotation_colors=annotation_colors,
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         annotation_colors=annotation_colors,
         annotation_row=annotation_row,
         annotation_col=annotation_col,
         main = "Zscore of mean exp genes for each cluster",scale="row")
dev.off()
