library(reshape2)
args=commandArgs(TRUE)
input=args[1]
level2=args[2]
level3=args[3]
fun=function(x) length(unique(unlist(strsplit(x,","))))
t=read.table(input,sep="\t",header=T,quote="!")

## level2
f=t[,c(1,2,6)]
colnames(f)=c("V1","V2","V3")
m = melt(f,id=c("V1","V2"))
d=dcast(m,V1+V2~variable,fun=fun)
colnames(d)=c("GO:level1","GO:level2","counts")
write.table(d,file=level2,sep="\t",col.names=T,row.names=F,quote=F)

## level3
f=t[,c(1,3,6)]
colnames(f)=c("V1","V2","V3")
m = melt(f,id=c("V1","V2"))
d=dcast(m,V1+V2~variable,fun=fun)
colnames(d)=c("GO:level1","GO:level3","counts")
write.table(d,file=level3,sep="\t",col.names=T,row.names=F,quote=F)

