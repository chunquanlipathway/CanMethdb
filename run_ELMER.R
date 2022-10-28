library(GenomicRanges)

# U.test: probe methylation value (front and back 20%) and near genes
# Load gene expression file
# Gene expression file: gene expression profile
GeneExp<-read.table("gene_expression.txt",header=T,sep="\t",row.names=1)
rownames(GeneExp)<-apply(as.matrix(rownames(GeneExp)),1,function(x){unlist(strsplit(x,"[.]"))}[1])
# Load main function
source("ELMER.R")

for(i in 1:22){
	# Load methylation profile splited by chromosome previously and gene expression profile
	Meth<-read.table(file=paste("methylation_chr",i,".txt",sep=""),header=T,row.names=1,sep="\t")
	mee <- list(meth=Meth,exp=GeneExp)
	rm(Meth)
	# U.test
	NearGenes<-read.table(file=paste("candidate.pair_chr",i,".bed",sep=""),header=T,sep="\t")
	Probe<-rownames(mee$meth)
	Gene<-rownames(mee$exp)
	nearGenes<-NearGenes[NearGenes$Target%in%Probe,]
	nearGenes<-nearGenes[nearGenes$GeneID%in%Gene,]
	rm(NearGenes)
	rm(Probe)
	rm(Gene)
	probes<-unique(nearGenes$Targe)##12782
	Probe.gene<-get.pair(mee,probes=probes,nearGenes=nearGenes,percentage=0.2,cores=4,portion = 0.3,save=FALSE)
	write.table(Probe.gene, file=paste("ELMER_chr",i,".bed",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
	print(i)

}











