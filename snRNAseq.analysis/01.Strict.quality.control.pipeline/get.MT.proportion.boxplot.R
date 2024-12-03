library(reshape2)
require("RColorBrewer")
setwd('/Users/zhangjiaman/Desktop/students.works/zjm.data/hic项目/34.Mouse.LDM.snRNAseq.for.HaoHuan/New.snRNAseq.analysis/01.snRNAseq.integrate.pipeline/raw.UMAP.for.MT.filter/')
num<-'2'
inp<-paste('whole.res.',num,'.median.MT.txt',sep='')
dat<-read.table(inp,head=F,sep='\t',stringsAsFactors = F)

pdf('cluster.median.MT.boxplot.pdf',w=2.3,h=4)

par(xaxs='i',yaxs='i',lend=2,mgp=c(1.4,0.6,0),tcl=-0.3,mar=c(5,5,4,2),bty='l',xpd=NA)
plot(1:10,type='n',axes=FALSE,xlab="",ylab="Proportion of MT reads (%)",xlim=c(0.7,1.5),ylim=c(0,13))

boxplot(dat$V2,at=1,lty=1,lwd=1, boxwex=0.7,add=TRUE,col='grey',staplewex=0.2,outline=FALSE,axes=F,xpd=NA)

out_list<-sort(boxplot(dat$V2,plot=F)$out,decreasing=T)
odat<-data.frame(x=rep(1,time=length(out_list)),out=out_list)

odat<-dat[dat$V2>=min(out_list),]
odat$x<-1

col_list<-rev(brewer.pal(8,'Accent')[1:6])[1:length(out_list)]
odat$col<-col_list
odat$title<-paste('MT-',odat$V1,sep='')

lines(odat$x,odat$V2,type='p',pch=16,col=odat$col,cex=0.6)

text(1.3,odat$V2,labels=odat$title,col=odat$col,cex=0.8)

axis(2,las=2)
box(bty='o')
dev.off()
