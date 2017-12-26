#!/bin/Rscript

args=commandArgs(T)
if(length(args) != 3) stop( "\tParameters count not correct! \n\n\tRscript plot.R <sample_name> <R1_file_prefix> <R2_file_prefix>\n\n")


sample_name=args[1]
R1_file_prefix=args[2]
R2_file_prefix=args[3]

title=sample_name

#*********************	GC_Content
out_GC=paste(sample_name,"-GCContent.pdf",sep="")
rt1=read.table(paste(R1_file_prefix,".GCContent",sep=""))
if(file.exists(paste(R2_file_prefix,".GCContent",sep=""))) {
	rt2=read.table(paste(R2_file_prefix,".GCContent",sep=""))
	rt=rt1+rt2
}else{
	rt=rt1
}

rt[[8]] = apply(rt[,2:6],1,sum)
pdf(out_GC,width=12,height=8)
color=rainbow(5)
length=length(rt[[1]])

plot(rt1[[1]],rt[[2]]*100/rt[[8]],type="l",lwd=2,cex=1.5,col=color[1],xlim=c(1,length),ylim=c(0,100),xlab="",ylab="",xaxt="n",main=title)
axis(1, at=seq(0,length,10),cex.axis=1)
lines(rt1[[1]],rt[[3]]*100/rt[[8]],type="l",lwd=2,cex=1.5,col=color[2])
lines(rt1[[1]],rt[[4]]*100/rt[[8]],type="l",lwd=2,cex=1.5,col=color[3])
lines(rt1[[1]],rt[[5]]*100/rt[[8]],type="l",lwd=2,cex=1.5,col=color[4])
lines(rt1[[1]],rt[[7]]*100/rt[[8]],type="l",lwd=2,cex=1.5,col=color[5])
title(,ylab="Percentage(%)",xlab="Position In Read(bp)",cex.main=2.0)
legend(length-30,100,lty = c(1,1,1,1,1),legend=c("A","C","G","T","GC"),cex=1.3,col=color)
dev.off()

#*********************	MeanQV Plot
out_MeanQV=paste(sample_name,"-MeanQVofRead.pdf",sep="")
pdf(out_MeanQV,width=12,height=8)

par(mfrow=c(2,1))
rt = read.table(paste(R1_file_prefix,".MeanQual",sep=""))
length=length(rt[[1]])
data=rt[,2]
pos =rt[,1]
names(data)=pos
barplot(data,col=rev(rainbow(length)),main=paste("Reads mean quality distribution of ",title,sep=""),xlab="Mean Quality Score of Read1",ylab="Reads Count") 
axis(side=1,at=pos,label=pos,xaxt="n")

if(file.exists(paste(R2_file_prefix,".MeanQual",sep=""))) {
rt = read.table(paste(R2_file_prefix,".MeanQual",sep=""))
length=length(rt[[1]])
data=rt[,2]
pos =rt[,1]
names(data)=pos
barplot(data,col=rainbow(length),xlab="Mean Quality Score of Read2",ylab="Reads Count") 
axis(side=1,at=pos,label=pos,xaxt="n")
}
dev.off()

#*********************	BaseQV BoxPlot
out_BaseQV=paste(sample_name,"-BaseQVofReads.pdf",sep="")
pdf(out_BaseQV,width=12,height=8)

par(mfrow=c(2,1))
rt = read.table(paste(R1_file_prefix,".QualDist",sep=""))
length=length(rt[,1])
data=t(rt)
summarydata = list(stats=data,n=rep(10,each=length))
bxp(summarydata,main=paste("Quality score across all bases of ",title,sep=""),xlab="Position In Read(bp) Read1",ylab="Quality Score")

if(file.exists(paste(R2_file_prefix,".QualDist",sep=""))) {
rt = read.table(paste(R2_file_prefix,".QualDist",sep=""))
length=length(rt[,1])
data=t(rt)
summarydata = list(stats=data,n=rep(10,each=length))
bxp(summarydata,xlab="Position In Read(bp) Read2",ylab="Quality Score")
}
dev.off()



