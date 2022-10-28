# Main function of finding near genes
NearGenes <- function (Target=NULL,Gene=NULL,geneNum=20,TRange=NULL){
	if(is.null(Gene) | is.null(Target)){
    stop ("Target and Genes should both be defined")
	}
	message(Target)
	if(is.null(TRange)){
		stop( "TRange must be defined")
	}else{
		regionInfo <- TRange[as.character(TRange$name) %in% Target]
	}
	GeneIDs <-c()
	Distances <- c()
	strand(Gene) <- "*"
	Gene <- Gene[as.character(seqnames(Gene)) %in% as.character(seqnames(regionInfo))]
	if(length(Gene)==0){
		warning(paste0(Target," don't have any nearby gene in the given gene list."))
		Final <- NA
	}else{
		Gene <- sort(Gene)
		index <- follow(regionInfo,Gene)
	
    # Find left side near genes
    Leftlimit <- geneNum/2
    Rightlimit <- geneNum/2
    n <- 1
    if(is.na(index)){
      index<- 0
      Leftlimit <- 0
      Left <- c()
    }else if(index==1){
      Left <- index
      Leftlimit <- length(Left)
    }else{
      Left <- index
      while(length(Left) < Leftlimit){
        if(!as.character(Gene$GENEID[index-n])%in%as.character(Gene$GENEID[Left]))
          Left <- c((index-n),Left) 
        if((index-n)==1){
          Leftlimit <- length(Left)
        } 
        n <- n+1
      }
    }
    
	# Find right side near genes
    Right <- c()
    n <- 1
    if(index==length(Gene) || 
         all(unique(Gene$GENEID[(index+1):length(Gene)]) %in% as.character(Gene$GENEID[index]))){
      Rightlimit <- length(Right)
    }else{
      while(length(Right) < Rightlimit){
        if(!as.character(Gene$GENEID[index+n])%in% as.character(Gene$GENEID[c(Right,Left)])) 
          Right <- c(Right,(index+n))
        
        if(index+n==length(Gene)){
          Rightlimit <- length(Right)
        } else{
          n <- n+1
        }     
      }
    }
   
    # If the near genes on the one side are insufficient, increase the number of near genes on the other side
    if(Rightlimit < geneNum/2){
      n <- 1
      if(Left[1]-n > 0){
        while((length(Left)+length(Right)) < geneNum){
          if(!as.character(Gene$GENEID[Left[1]-n])%in%as.character(Gene$GENEID[c(Left,Right)])) 
            Left <- c((Left[1]-n),Left) 
          n <- n+1
        }
      }  
    }
    
    if(Leftlimit < geneNum/2){
      n <- 1
      m <- length(Right)
      if(Right[m]+n < length(Gene)+1)
        while((length(Left)+length(Right)) < geneNum){
          if(!as.character(Gene$GENEID[Right[m]+n])%in%as.character(Gene$GENEID[c(Left,Right)])) 
            Right <- c(Right,(Right[m]+n)) 
          n <- n+1
        } 
    }
    Whole <- c(Left,Right)
    GeneIDs <- Gene$GENEID[Whole]
    Symbols <- Gene$SYMBOL[Whole]
    Distances <-  suppressWarnings(distance(Gene[Whole],regionInfo))
    if(Rightlimit < 1){
      Sides <- paste0("L",length(Left):1)
    }else if( Leftlimit < 1){
      Sides <- paste0("R",1:length(Right))
    }else{
      Sides <- c(paste0("L",length(Left):1),paste0("R",1:length(Right)))
    }
    Final <- data.frame(Target=rep(Target,length(GeneIDs)),GeneID=GeneIDs,
                        Symbol=Symbols,Distance=Distances, Side=Sides, 
                        stringsAsFactors = FALSE)
  }
  return(Final)
}