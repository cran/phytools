# function to perform phylogenetic principal components analysis
# multiple morphological traits in Y
# also can use lambda transformation in which lambda is optimized by ML
# written by Liam Revell 2010/2011, ref. Revell (2009; Evolution)

phyl.pca<-function(tree,Y,method="BM",mode="cov"){
	# check tree
	if(class(tree)!="phylo") stop("tree must be an object of class 'phylo.'")
	# check mode
	if(length(strsplit(mode,split="")[[1]])<=2){ 
		message(paste("mode = '",mode,"' not a valid option; setting mode = 'cov'",sep=""))
		mode="cov"
	}
	if(all(strsplit(mode,split="")[[1]]==strsplit("correlation",split="")[[1]][1:length(strsplit(mode,split="")[[1]])]))
		mode="corr"
	else if(all(strsplit(mode,split="")[[1]]==strsplit("covariance",split="")[[1]][1:length(strsplit(mode,split="")[[1]])]))
		mode="cov"
	else {
		message(paste("mode = '",mode,"' not a valid option; setting mode = 'cov'",sep=""))
		mode="cov"
	}
	# preliminaries
	n<-nrow(Y); m<-ncol(Y)
	# check and sort data
	if(n>length(tree$tip)) stop("number of rows in Y cannot be greater than number of taxa in your tree")
	Y<-as.matrix(Y)
	if(is.null(rownames(Y))){
		if(nrow(Y)==n){ 
			print("Y has no names. function will assume that the row order of Y matches tree$tip.label")
			rownames(Y)<-tree$tip.label
		} else stop("Y has no names and does not have the same number of rows as tips in tree")
	} else if(length(setdiff(rownames(Y),tree$tip.label))!=0) 
		stop("Y has rownames, but some rownames of Y not found in tree")
	# analyze
	C<-vcv.phylo(tree)[rownames(Y),rownames(Y)]
	if(method=="BM"){ 
		temp<-phyl.vcv(Y,C,1.0)
		V<-temp$R; a<-t(temp$alpha); C<-temp$C
	} else if(method=="lambda"){
		temp<-optimize(f=likMlambda,interval=c(0,maxLambda(tree)),X=Y,C=C,maximum=TRUE)
		lambda<-temp$maximum; logL<-as.numeric(temp$objective)
		temp<-phyl.vcv(Y,C,lambda)
		V<-temp$R; a<-t(temp$alpha); C<-temp$C
	}
	invC<-solve(C)
	# if correlation matrix
	if(mode=="corr"){
		Y=Y/matrix(rep(sqrt(diag(V)),n),n,m,byrow=T) # standardize Y
		V=V/(sqrt(diag(V))%*%t(sqrt(diag(V)))) # change V to correlation matrix
		a<-matrix(colSums(invC%*%Y)/sum(invC),m,1) # recalculate a
	}
	es=eigen(V) # eigenanalyze
	result<-list(); result$Eval<-diag(es$values); result$Evec<-es$vectors
	dimnames(result$Eval)<-list(paste("PC",1:ncol(Y),sep=""),paste("PC",1:ncol(Y),sep=""))
	dimnames(result$Evec)<-list(colnames(Y),paste("PC",1:ncol(Y),sep=""))
	A<-matrix(rep(a,n),n,m,byrow=T)
	result$S<-(Y-A)%*%result$Evec # compute scores in the species space
	Ccv<-t(Y-A)%*%invC%*%result$S/(n-1) # compute cross covariance matrix and loadings
	result$L<-matrix(,m,m,dimnames=list(colnames(Y),paste("PC",1:ncol(Y),sep="")))
	for(i in 1:m) for(j in 1:m) result$L[i,j]<-Ccv[i,j]/sqrt(V[i,i]*result$Eval[j,j])
	if(method=="lambda"){ 
		result$lambda<-lambda
		result$logL.lambda<-logL
	}
	# return result
	return(result)
}


