library(GenomicRanges)
# Load main function
source('NearGenes.R')

# Split probe file by chromosome
# Probe file: Chr, Start position, Probe name 
ProbeLoc<-read.table("probe.bed",header=F,sep=" ")
chr_i<-substr(x=ProbeLoc[,1],start=4,stop=5)
probe<-cbind(ProbeLoc,chr_i)
colnames(probe)<-c("chr","Start","name","chr_i")
for(i in 1:22){
	probe_chri<-probe[which(probe[,4]==i),]
	write.table(probe_chri[,1:3],file=paste("probe_chr",i,".bed",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
}

# Split gene file by chromosome
# Gene file: Chr, Start position, End position, Gene id, Gene SYMBOL id
GeneLoc<-read.table("gene.txt",header=F,sep="\t")
chr_i<-substr(x=GeneLoc[,1],start=4,stop=5)
gene<-cbind(GeneLoc,chr_i)
colnames(gene)<-c("chr","start","end","id","SYMBOL","chr_i")
for(i in 1:22){
	gene_chri<-gene[which(gene[,6]==i),]
	write.table(gene_chri[,1:5],file=paste("gene_chr",i,".bed",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
}

# Parallel preparation
cores<-parallel::detectCores()
if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
		if(!is.null(cores)){
			if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
			cl <- snow::makeCluster(cores,type = "SOCK")
		}
	}

# Run in parallel by chromosome
for(i in 1:22){
	probe_chri<-read.table(file=paste("probe_chr",i,".bed",sep=""),sep="\t",header=T)
	TRange<- GRanges(seqnames=probe_chri$chr,ranges=IRanges(probe_chri$Start,width=1),name=probe_chri$name)
	gene_chri<-read.table(file=paste("gene_chr",i,".bed",sep=""),sep="\t",header=T)
	Gene<-GRanges(seqnames=gene_chri$chr,ranges=IRanges(gene_chri$start,gene_chri$end,width=0),tx_name=gene_chri$SYMBOL,GENEID=gene_chri$id,SYMBOL=gene_chri$SYMBOL)
	NearGene_i<- parallel::parSapplyLB(cl,as.character(TRange$name),NearGenes,
			                                     geneNum=20,Gene=Gene,TRange=TRange,
			                                     simplify=FALSE)
	# If you cannot run in parallel, run sapply									
	# NearGene_i<- sapply(as.character(TRange$name),NearGenes,geneNum=20,
	#		                     Gene=Gene,TRange=TRange,simplify=FALSE)
	NearGene<-do.call(rbind,NearGene_i)
	write.table(NearGene,file=paste("candidate.pair_chr",i,".bed",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
	print(i)
}
