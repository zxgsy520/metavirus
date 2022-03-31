#!/usr/bin/env Rscript

options(bitmapType='cairo')

args=commandArgs(TRUE)
input=args[1]
output=args[2]
options(stringsAsFactors=F)
t=read.table(input,header=T,sep="\t",quote="!",comment.char="") #GO.annotion.xls
#N = 100

col = c("lightseagreen","lightcoral","mediumpurple")
#t[,3] = log(t[,3],2)
sum=sum(t[,3])
t[,3] = t[,3]/sum*100
t[,3] = log(t[,3],10)
max = max(t[,3])
min = min(t[,3])
#t=t[t[,4]>=min,-5]
mf = "molecular_function"
bp = "biological_process"
cc = "cellular_component"
mbc=sapply(strsplit(t[,1],"\\(|\\)"),'[',2)

MF = t[mbc==mf,]
n1 = nrow(MF)
BP = t[mbc==bp,]
n2 = nrow(BP)
CC = t[mbc==cc,]
n3 = nrow(CC)
y = c(MF[,3],BP[,3],CC[,3])
labels = c(MF[,2],BP[,2],CC[,2])
color = rep(col,c(n1,n2,n3))
pdf(output,width=12)
par(mar=c(15,8,4,3))
d=ifelse(sum(y< -2)>5,2,1)
y0=y+d
max=max+d
min=min+d
y0=ifelse(y0>0,y0,0)
par(las=1,mgp=c(3,0.5,0),cex.axis=0.6)
bar=barplot(y0,col=color,border=color,space=1,ylab="Percent of Genes",ylim=c(-1*d,2)+d,width=1.5,yaxt="n")
par(xpd=T,srt=45)
text(bar,-0.2,labels=labels,col=color,cex=0.5,adj=1)
par(srt=0)
axis(2,(-1*d):2+d,10^((-1*d):2))
axis(4,(-1*d):2+d,round(10^((-1*d):2)*sum/100,0))
text((0+bar[n1])/2,max*1.1,mf,cex=0.7,col=col[1])
text((bar[n1]+bar[n1]+bar[n2])/2,max*1.1,bp,cex=0.7,col=col[2])
text((bar[n1]+bar[n2]+bar[n1]+bar[n2]+bar[n3])/2,max*1.1,cc,cex=0.7,col=col[3])
par(las=3)
mtext("Number Of Genes",side=4,line=1.5)
dev.off()
