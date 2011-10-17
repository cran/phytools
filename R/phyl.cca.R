# function does phylogenetic canonical correlation analysis (Revell & Harrison 2008)
# written by Liam Revell 2011

phyl.cca<-function(tree,X,Y,lambda=1.0,fixed=TRUE){
	# misc
	n<-length(tree$tip)
	mX<-ncol(X)
	mY<-ncol(Y)
	# compute C
	C<-vcv.phylo(tree)
	# reorder Y & C by X
	C<-C[rownames(X),rownames(X)]
	Y<-Y[rownames(X),]
	# set or optimize lambda
	if(fixed) C<-lambda.transform(lambda,C)
	else {
		temp<-optimize(f=likelihood,interval=c(0,1),X=cbind(X,Y),C=C,maximum=TRUE)
		lambda<-temp$maximum
		C<-lambda.transform(lambda,C)
	}
	# invert C
	invC<-solve(C)
	# compute means
	aX<-colSums(invC%*%X)/sum(invC)
	aY<-colSums(invC%*%Y)/sum(invC)
	# compute cov & cross-cov matrices
	one<-as.matrix(rep(1,n))
	SigXX<-t(X-one%*%aX)%*%invC%*%(X-one%*%aX)/(n-1)
	SigXY<-t(X-one%*%aX)%*%invC%*%(Y-one%*%aY)/(n-1)
	SigYX<-t(SigXY)
	SigYY<-t(Y-one%*%aY)%*%invC%*%(Y-one%*%aY)/(n-1)
	# compute canonical coefficients
	A<-eigen(solve(SigXX)%*%SigXY%*%solve(SigYY)%*%SigYX)
	B<-eigen(solve(SigYY)%*%SigYX%*%solve(SigXX)%*%SigXY)
	# compute canonical variables, rescale
	U<-X%*%A$vectors[,1:min(mX,mY)]
	aU<-colSums(invC%*%U)/sum(invC)
	vcvU<-t(U-one%*%aU)%*%invC%*%(U-one%*%aU)/(n-1)
	U<-(U-one%*%aU)%*%diag(sqrt(diag(1/vcvU)))
	V<-Y%*%B$vectors[,1:min(mX,mY)]
	aV<-colSums(invC%*%V)/sum(invC)
	vcvV<-t(V-one%*%aV)%*%invC%*%(V-one%*%aV)/(n-1)
	V<-(V-one%*%aV)%*%diag(sqrt(diag(1/vcvV)))
	# compute canonical correlations
	aU<-colSums(invC%*%U)/sum(invC)
	aV<-colSums(invC%*%V)/sum(invC)
	Ccv<-round(t(cbind(U,V))%*%invC%*%cbind(U,V)/(n-1),10)
	ccs<-diag(Ccv[1:min(mX,mY),(1+min(mX,mY)):(2*min(mX,mY))])
	pos<-2*(as.numeric(ccs>0)-0.5)
	ccs<-ccs*pos
	# reorient variables, reorient & rescale coefficents
	U<-U*one%*%pos
	xcoef<-A$vectors[,1:min(mX,mY)]*matrix(1,mX,1)%*%pos%*%diag(sqrt(diag(1/vcvU)))
	ycoef<-B$vectors[,1:min(mX,mY)]%*%diag(sqrt(diag(1/vcvV)))
	# conduct hypothesis tests
	W_lh<-rep(1,min(mX,mY))
	chiSq<-vector()
	df<-vector()
	for(i in 1:min(mX,mY)){ 
		for(j in i:min(mX,mY)) W_lh[i]<-W_lh[i]*(1-ccs[j]^2)
		chiSq[i]<--((n-1)-(mX+mY+1)/2)*log(W_lh[i])
		df[i]<-(mX+1-i)*(mY+1-i)
	}
	pvalues<-pchisq(chiSq,df=df,lower.tail=F)
	# add row & column names
	if(!is.null(colnames(X))) rownames(xcoef)<-colnames(X)
	if(!is.null(colnames(Y))) rownames(ycoef)<-colnames(Y)
	temp<-vector()
	for(i in 1:min(mX,mY)) temp[i]<-paste("CA",i,sep="")
	colnames(xcoef)<-temp
	colnames(ycoef)<-temp
	colnames(U)<-temp
	colnames(V)<-temp
	# return as list
	return(list(cor=ccs,xcoef=xcoef,ycoef=ycoef,xscores=U,yscores=V,lambda=lambda,chisq=chiSq,p=pvalues))
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
	
	
