library(DMwR)
library(limma)

for(i in 1:22){
    
	#load CpG-target gene pairs computed by the two methods  
	cor_result <- read.table(file=paste("CpG-gene pairs_chr",i,".txt",sep=""),sep="\t",header=T)
  	#
	if(dim(cor_result)[1]==0) {
		next
	}
		else{
			probe <- unique(cor_result$Keys)
			exprSet <- read.table(file=paste("methylation_chr",i,".txt",sep=""),sep="\t",header=T,row.names= 1)
			exprSet <- exprSet[which(rownames(exprSet)%in%probe),]
			# Group by disease and normal samples
			samplenames <- as.matrix(colnames(exprSet))
			group_list2 <- as.matrix(as.numeric(substr(samplenames,start=14,stop=15)))
			group_list <- as.matrix(apply(group_list2,1,function(x){if(x <= 9) group_list <- "case" else group_list<- "control"}))
			design <-  model.matrix(~0+factor(group_list))
			colnames(design)=as.character(levels(factor(group_list)))
			rownames(design)=colnames(exprSet)
			# Establish regression model and group the models to find differences
			contrast.matrix <- makeContrasts(paste0(colnames(design),collapse = "-"),levels = design)
			fit3 <-  lmFit(exprSet,design)
			fit2 <-  contrasts.fit(fit3, contrast.matrix)
			fit <-  eBayes(fit2)
			tempOutput <-  topTable(fit, coef=1, n=Inf)
			tempOutput <- tempOutput[,4:5]
			minus <- as.data.frame(rowMeans(exprSet[,which(group_list=="case")])-rowMeans(exprSet[,which(group_list=="control")]))
			colnames(minus) <- "minus"
			minus_re <- as.data.frame(minus[match(row.names(tempOutput),row.names(minus)),])
			colnames(minus_re) <- "minus"
			output <- cbind(tempOutput,minus_re)
			minus_re[minus_re>0]="hyper"
			minus_re[minus_re<0]="hypo"
			out.with.lable<-cbind(output ,minus_re)
			limma_0.1<-out.with.lable[which(abs(out.with.lable$minus)>=0.1),]
			colnames(limma_0.1)<-c("limma.P.Value","limma.adj.P.Val","limma.minus","lable")
			write.table(limma_0.1,file=paste("limma_chr",i,".txt",sep=""),sep="\t",row.names=T,col.names=T,quote=F)
		}

}