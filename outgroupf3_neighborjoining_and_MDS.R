## Outgroupf3results.par.log is the outgroupf3 results file (pruned to just the row of results from ADMIXTOOLS).
AllIndividuals = read.table("Outgroupf3results.par.log",header=F)
AllIndividuals$V1 = as.character(AllIndividuals$V1)
AllIndividuals$V3 = as.character(AllIndividuals$V3)
Allpops=unique(AllIndividuals$V2)
Allpops=as.character(Allpops)
## Add last pop to list
test=as.character(AllIndividuals[nrow(AllIndividuals),3])
Allpops=c(Allpops,test)
Table=data.frame(1:length(Allpops))
for(i in 1:length(Allpops)){
Table[,i]=1
}
c=1
for(i in 1:(length(Allpops))){
for(j in (i+1):(length(Allpops))){
Table[i,j]=AllIndividuals[c,5]
c=c+1
}}
Table[,ncol(Table)]=NULL
c=1
for(i in 1:(length(Allpops))){
for(j in (i+1):(length(Allpops))){
Table[j,i]=AllIndividuals[c,5]
c=c+1
}}
Table <- Table[-nrow(Table),]
Table[nrow(Table),nrow(Table)]=1

### make middle columns 0 for neighbor joining tree.
for(i in 1:(nrow(Table))){
Table[i,i]=0
}
write.table(Table,file="outgroupf3_forneighborjoiningtree.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

####### Make a distance matrix for MDS plot.
Table=1-Table
### make middle columns 0.
for(i in 1:(nrow(Table))){
Table[i,i]=0
}
write.table(Table,file="outgroupf3_forMDS.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)


## Prepare .phyl file for Itol plotting.
AllIndividuals = read.table("Outgroupf3results.par.log",header=F)
AllIndividuals$V1 = as.character(AllIndividuals$V1)
AllIndividuals$V3 = as.character(AllIndividuals$V3)
Allpops=unique(AllIndividuals$V2)
Allpops=as.character(Allpops)
## Add last pop to list
test=as.character(AllIndividuals[nrow(AllIndividuals),3])
Allpops=c(Allpops,test)
Table=read.table("outgroupf3_forneighborjoiningtree.txt",header=F)
Table=1/Table
### make middle columns 0.
for(i in 1:(nrow(Table))){
Table[i,i]=0
}
Table2=data.frame(1:ncol(Table))
odd=sapply(Allpops,function(x) substr(x,1,10))
Table2[,1]=odd
Table3=cbind(Table2,Table)
write.table(Table3,file="outgroupf3_forneighborjoiningtree.phyl",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)


##### Plotting MDS results.
dist.au <- read.table("outgroupf3_forMDS.txt")
fit <- cmdscale(dist.au, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
pdf("outgroupf3_MDSplot.pdf",width=8,height=6.5, useDingbats=FALSE)
city.names= read.table("NamesofGroups.txt",header=F)
city.names$V1=as.character(city.names$V1)
city.names$V2=as.character(city.names$V2)
plot(x,y, pch = 19, ylim=c(-0.23,0.22), xlim=c(-0.2,0.2),col=city.names$V2)
text(x,y, pos = 4, labels = city.names$V1, cex=0.5)
dev.off()
