# Main function and subfunction
# Function of U.test in each probe
Stat.nonpara <- function(Probe,NearGenes,K,Top=NULL,Meths=Meths,Exps=Exps){
	if(! length(Probe)==1) {stop("Number of Probe should be 1")}
	# Take out the near gene corresponding to the probe
	Gene<-NearGenes[which(NearGenes[,1]%in%Probe),]$GeneID
	# Take out the corresponding expression values of near genes
	Exp <- Exps[as.character(Gene),]
	Meth <- Meths[Probe,]
	# Binary
	Meth_B <- mean(Binary(as.numeric(Meth),Break=K),na.rm = TRUE)
	if( Meth_B >0.95 | Meth_B < 0.05 ){
		test.p <- NA
	}else{
	# Take out the 20% sample expression value methy and unmethy
    unmethy <- order(Meth)[1:round(length(Meth)*Top)]
    methy <- order(Meth,decreasing=TRUE)[1:round(length(Meth)*Top)] 
    Fa <- factor(rep(NA,length(Meth)),levels=c(-1,1))
    Fa[unmethy] <- -1
    Fa[methy] <- 1
    if(!is.vector(Exp)){
		Exp <- as.matrix(Exp[,!is.na(Fa)])
		Fa <- Fa[!is.na(Fa)]
		test.p <- unlist(lapply(splitmatrix(Exp),function(x,Factor){
													wilcox.test(x[Factor %in% -1],x[Factor %in% 1],
													alternative = "greater",
													exact=FALSE)$p.value
													},
						Factor=Fa))
    }else{
		Exp <- as.matrix(Exp[!is.na(Fa)])
		Fa <- Fa[!is.na(Fa)]
		test.p <- wilcox.test(Exp[Fa %in% -1],Exp[Fa %in% 1],
							  alternative = "greater",
                              exact=FALSE)$p.value
		}
	}
 
	if(length(Gene)==1){
		out <- data.frame(Probe=rep(Probe,length(Gene)),
						  GeneID=Gene,Symbol=NearGenes[which(NearGenes[,1]%in%Probe),]$Symbol, 
                          Distance=NearGenes[which(NearGenes[,1]%in%Probe),]$Distance, 
                          Sides=NearGenes[which(NearGenes[,1]%in%Probe),]$Side,
                          Raw.p=test.p, 
                          stringsAsFactors = FALSE)
	}else{
		out <- data.frame(Probe=rep(Probe,length(Gene)),
                          GeneID=Gene,Symbol=NearGenes[which(NearGenes[,1]%in%Probe),]$Symbol, 
                          Distance=NearGenes[which(NearGenes[,1]%in%Probe),]$Distance, 
                          Sides=NearGenes[which(NearGenes[,1]%in%Probe),]$Side,
                          Raw.p=test.p[match(Gene, names(test.p))], 
                          stringsAsFactors = FALSE)
	}
	print(Probe)
	return(out)
}

Binary <- function(x,Break=0.3,Break2=NULL){
	if(!is.numeric(x)) stop("x need to be numeric") 
	change <- x
	if(is.null(Break2)){
		change[x > Break] <- 1
		change[x < Break | x== Break] <- 0
	}else{
		change[x < Break | x== Break] <- 0
		change[x> Break & x < Break2] <- NA
		change[x > Break2 | x== Break2] <-1 
	} 
	return(change)    
}

splitmatrix <- function(x,by="row"){
	if(by %in% "row"){
		out <- split(x, rownames(x))
	}else if (by %in% "col"){
		out <- split(x, colnames(x))
	}
	return(out)
}

count<-function(x){
	re<-c(length(which(x>Pp)))
}

# Function of disturbance U.test in each probe
Stat.nonpara.permu <- function(Probe,Gene,Top=0.2,Meths=Meths,Exps=Exps,permu.dir=NULL){
	Exp <- Exps[Gene,]
	if(is.vector(Meths)){
		Meth <- Meths
	}else{
		Meth <- Meths[Probe,]
	}
	unmethy <- order(Meth)[1:round(length(Meth)*Top)] 
	methy <- order(Meth,decreasing=TRUE)[1:round(length(Meth)*Top)] 
	Fa <- factor(rep(NA,length(Meth)),levels=c(-1,1))
	Fa[unmethy] <- -1
	Fa[methy] <- 1
	Exp <- Exp[!is.na(Fa)]
	Fa <- Fa[!is.na(Fa)]
	test.p <-wilcox.test(as.numeric(Exp[Fa %in% -1]),as.numeric(Exp[Fa %in% 1]),alternative = "greater",exact=FALSE)$p.value
	out <- data.frame(Gene,Probe,test.p, 
					  stringsAsFactors = FALSE) 
	return(out)
}

# Random disturbance
get.permu <- function(geneID, mee, percentage=0.2, portion=0.3,
                      permu.size=permu.size){
	# Set seed
	set.seed(200)
	# get usable probes
	rm.probes <- Probe.gene[Probe.gene$GeneID%in%geneID,1]
	binary.m <- rowMeans(Binary(as.matrix(mee$meth),portion),na.rm = TRUE)
	usable.probes <- names(binary.m[binary.m <0.95 & binary.m > 0.05 & !is.na(binary.m)])
	usable.probes <- usable.probes[!usable.probes %in% rm.probes]
	if(length(usable.probes) < permu.size) 
		stop(sprintf("There is no enough usable probes to perform %s time permutation, 
					set a smaller permu.size.",permu.size))
	if(!is.numeric(permu.size)) permu.size <- length(usable.probes) 
	probes.permu <- as.matrix(sample(usable.probes, size = permu.size, replace = FALSE))
	permu.meth<-mee$meth[match(probes.permu,rownames(mee$meth)),]
	permu<-do.call(rbind,sapply(probes.permu,Stat.nonpara.permu,Meths=permu.meth,
								Gene=geneID,
								Top=0.2,Exps=mee$exp,
								simplify=FALSE))
	Pr <- Probe.gene[Probe.gene$GeneID%in%geneID,6]
	Pp <- permu$test.p
	num<-apply(as.matrix(Pr),1,count)
	Pe <- (num+1)/permu.size
	permu.Pe <- as.matrix(cbind(Probe.gene[Probe.gene$GeneID%in%geneID,],Pe))
	print(geneID)
	return(permu.Pe)
}

# Main function
get.pair <- function(mee,probes,nearGenes,percentage=0.2,permu.size=NULL,
                     permu.dir=NULL,dir.out="./",cores=NULL,
                     portion = 0.3,save=TRUE){
	# check data
	if(!all(probes %in% rownames(mee$meth))) 
		stop("Probes option should be subset of rownames of methylation matrix.")
	if(is.character(nearGenes)){
		newenv <- new.env()
		load(nearGenes, envir=newenv)
		# The data is in the one and only variable
		nearGenes <- get(ls(newenv)[1],envir=newenv) 
	}else if(!is.list(nearGenes)){
		stop("nearGene option must be a list containing output of GetNearGenes function 
			or path of rda file containing output of GetNearGenes function.")
	}
	#get raw pvalue
	##I need to modify that if there is all NA. stop the process.
	if(requireNamespace("parallel", quietly=TRUE) && requireNamespace("snow", quietly=TRUE)) {
		if(!is.null(cores)){
			if(cores > parallel::detectCores()) cores <- parallel::detectCores()/2
			cl <- snow::makeCluster(cores,type = "SOCK")
			Probe.gene<-parallel::parSapplyLB(cl,as.matrix(probes),Stat.nonpara,Meths= mee$meth[rownames(mee$meth)%in%as.matrix(probes),], 
											  NearGenes=nearGenes,K=portion,Top=percentage,
											  Exps=mee$exp,simplify = FALSE)
			parallel::stopCluster(cl)
		}else{
			Probe.gene<-sapply(as.matrix(probes),Stat.nonpara,Meths=mee$meth[rownames(mee$meth)%in%as.matrix(probes),],
	                           NearGenes=nearGenes,K=portion,Top=percentage,Exps=mee$exp,
	                           simplify = FALSE)
		}
	} else {
		Probe.gene<-sapply(as.matrix(probes),Stat.nonpara,Meths=mee$meth[rownames(mee$meth)%in%as.matrix(probes),],
						   NearGenes=nearGenes,K=portion,Top=percentage,Exps=mee$exp,
                           simplify = FALSE)
	}
  
	Probe.gene <- do.call(rbind,Probe.gene)
	Probe.gene <- Probe.gene[!is.na(Probe.gene$Raw.p),]
	GeneID <- unique(Probe.gene[!is.na(Probe.gene$Raw.p),"GeneID"])
	if(save){
		write.table(Probe.gene,file=sprintf("getPair.pairs.significant.txt",dir.out), 
					row.names=FALSE,col.names=T,quote=F,sep="\t")
		}
	return(Probe.gene)
}

