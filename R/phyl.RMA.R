# this function computes a phylogenetic reduced major axis (RMA) regression
# written by Liam Revell 2010,2011,2012

phyl.RMA<-function(x,y,tree,method="BM",lambda=NULL,fixed=FALSE,h0=1.0){

	x<-x[tree$tip.label]; y<-y[tree$tip.label]

	# bind the x & y into columns
	X<-cbind(x,y)

	# function to compute C with lambda
	lambda.transform<-function(C,lambda){
		if(lambda==1) return(C)
		else {
			V<-diag(diag(C))
			C<-C-V
			C.lambda<-(V+lambda*C)
			return(C.lambda)
		}
	}

	# function to compute phylogenetic VCV using joint Pagel's lambda
	phyl.vcv<-function(X,phy,lambda){
		C<-vcv(phy)
		C<-lambda.transform(C,lambda)
		one<-matrix(1,nrow(C),1)
		a<-solve(t(one)%*%solve(C)%*%one)%*%(t(one)%*%solve(C)%*%X)
		V<-t(X-one%*%a)%*%solve(C)%*%(X-one%*%a)/nrow(C)
		return(list(R=V,alpha=a))
	}

	# likelihood function
	likelihood<-function(theta,X,phy){
		C<-vcv(phy)
		C<-lambda.transform(C,theta)
		# compute R, conditioned on lambda
		temp<-phyl.vcv(X,phy,theta);
		R<-temp$R; a<-temp$alpha; rm(temp)
		# prep
		n<-nrow(X); m<-ncol(X); D<-matrix(0,n*m,m)
		for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
		one<-matrix(1,n,1)
		y<-as.matrix(as.vector(X))
		# compute the log-likelihood	
		logL<--t(y-D%*%t(a))%*%solve(kronecker(R,C))%*%(y-D%*%t(a))/2-n*m*log(2*pi)/2-determinant(kronecker(R,C),logarithm=TRUE)$modulus[1]/2
		return(logL)
	}

	if(method=="lambda")
		if(fixed==FALSE)
			result<-optimize(f=likelihood,interval=c(0,1),X=X,phy=tree,maximum=TRUE)
		else
			result<-list(objective=likelihood(lambda,X,tree),maximum=lambda)
	else if(method=="BM")
		result<-list(objective=likelihood(1.0,X,tree),maximum=1.0)
	else
		stop("do not recognize method")	

	est.lambda<-result$maximum # estimated lambda

	C<-vcv(tree)
	C<-lambda.transform(C,est.lambda)

	temp<-phyl.vcv(X,phy=tree,lambda=est.lambda)

	beta1<-sqrt(temp$R[2,2]/temp$R[1,1])

	beta0<-temp$a[2]-beta1*temp$a[1]

	r<-y-(beta0+beta1*x)

	r2<-temp$R[1,2]^2/(temp$R[1,1]*temp$R[2,2])

	T<-abs(log(beta1)-log(h0))/sqrt((1-r2)/(length(tree$tip)-2))

	df<-2+(length(tree$tip)-2)/(1+0.5*r2)

	P<-2*pt(T,df=df,lower.tail=FALSE)

	test<-c(r2,T,df,P); names(test)<-c("r2","T","df","P")

	return(list(RMA.beta=c(beta0,beta1),V=temp$R,lambda=est.lambda,logL=as.numeric(result$objective),test=test,resid=as.matrix(r)))

}
