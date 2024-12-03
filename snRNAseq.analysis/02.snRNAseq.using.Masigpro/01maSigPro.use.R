#export R_LIBS=/Lustre01/tangqianzi/software/Rlibs35another/:$R_LIBS
#/Lustre02/00.software/R-3.5.0/bin/R

#================ run TPM start =======================================

library(maSigPro)
setwd('/Users/zhangjiaman/Desktop/students.works/zjm.data/hic项目/34.Mouse.LDM.snRNAseq.for.HaoHuan/New.snRNAseq.analysis/02.snRNAseq.downstram.analysis/Rename/07.pseudoBulk.for.each.cell.type/01.Masigpro/Immune/')
celltype<-'Immune'
data<-read.table(file=paste("../Matrix.data/",celltype,".DEGs.mean.UMI.matrix.txt",sep=''),sep="\t",header=T,check.names = F)
dim(data)
data<-na.omit(data)
dim(data)
head(data)
rownames(data)<-data$gid
data<-data[,-c(1)]
#data_design<-read.table(file="group_factor_file.txt",sep="\t",header=T,row.names=1)
data_design<-data.frame(Time=c(1,1,2,2,3,3),
                        Replicates=c(1,1,2,2,3,3),Group=c(1,1,2,2,3,3))
rownames(data_design)<-colnames(data)
#data_design<-data.frame(name=colnames(data),sample=c(rep(1,time=10),rep(2,time=6),rep(3,time=4),rep(4,time=6),rep(5,time=12)))
d = make.design.matrix(as.matrix(data_design), degree = 2)

fit = p.vector(data, d, Q=0.05, MT.adjust = "BH",counts = FALSE)
tstep = T.fit(fit,step.method="backward")

sigs = get.siggenes(tstep, rsq = 0.6, vars = "groups")

pdf("Figure 1.pdf")
cl = see.genes(sigs$sig.genes$Group,newX11 = F,alfa = 0.05, k = 4)## 对差异结果进行聚类分析
dev.off()
cluster = as.data.frame(cl$cut)  ## 提取聚类信息
head(cluster)

write.table(cluster,'Cluster.info.txt',sep='\t',col.names = F,quote=F)

