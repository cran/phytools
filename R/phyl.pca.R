# function to perform phylogenetic principal components analysis
# multiple morphological traits in Y
# also can use lambda transformation in which lambda is optimized by ML
# written by Liam Revell 2010/2011, ref. Revell (2009; Evolution)

phyl.pca<-function(tree,Y,method="BM",mode="cov"){
	# check 'ape'
	if(!require(ape)) stop("function requires 'ape' package. please install.")
	# check tree
	if(class(tree)!="phylo") stop("tree must be an object of class 'phylo.'")

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
		temp<-optimize(f=likelihood,interval=c(0,1),X=Y,C=C,maximum=TRUE)
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

# function to compute phylogenetic VCV using joint Pagel's lambda
# written by Liam Revell 2011
phyl.vcv<-function(X,C,lambda){
	C<-lambda.transform(lambda,C)
	invC<-solve(C)
	a<-matrix(colSums(invC%*%X)/sum(invC),ncol(X),1)
	A<-matrix(rep(a,nrow(X)),nrow(X),ncol(X),byrow=T)
	V<-t(X-A)%*%invC%*%(X-A)/(nrow(C)-1)
	return(list(C=C,R=V,alpha=a))
}

# lambda transformation of C
# written by Liam Revell 2011
lambda.transform<-function(lambda,C){
	if(lambda==1) return(C)
	else {
		V<-diag(diag(C))
		C<-C-V
		C.lambda<-(V+lambda*C)
		return(C.lambda)
	}
}

# likelihood function
# written by Liam Revell 2011
likelihood<-function(theta,X,C){
	# compute R, conditioned on lambda
	temp<-phyl.vcv(X,C,theta);
	C<-temp$C; R<-temp$R; a<-temp$alpha
	# prep
	n<-nrow(X); m<-ncol(X); D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	y<-as.matrix(as.vector(X))
	# compute the log-likelihood
	kronRC<-kronecker(R,C)
	logL<--t(y-D%*%a)%*%solve(kronRC,y-D%*%a)/2-n*m*log(2*pi)/2-log(det(kronRC))/2
	return(logL)
}


