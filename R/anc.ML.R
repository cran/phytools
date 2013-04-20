# lightweight version of ace(...,method="ML") for continuous traits
# written by Liam J. Revell 2011-2013

anc.ML<-function(tree,x,maxit=2000){

	# check to make sure that C will be non-singular
	if(any(tree$edge.length<=(10*.Machine$double.eps)))
		stop("some branch lengths are 0 or nearly zero")

	# function returns the log-likelihood
	likelihood<-function(par,C,invC,detC,x){
		sig2<-par[1]; a<-par[2]; y<-par[3:length(par)]
		z<-c(x,y)-a
		logLik<--z%*%invC%*%z/(2*sig2)-nrow(C)*log(2*pi)/2-nrow(C)*log(sig2)/2-detC/2
		return(-logLik)
	}

	# compute C
	C<-vcvPhylo(tree)
	invC<-solve(C)
	detC<-determinant(C,logarithm=TRUE)$modulus[1]
	x[rownames(C)[1:length(tree$tip)]]->x

	# assign starting values
	zz<-fastAnc(tree,x)
	y<-zz[2:length(zz)]
	a<-zz[1]
	sig2<-(c(x,y)-a)%*%invC%*%(c(x,y)-a)/nrow(C)

	# optimize
	res<-optim(c(sig2,a,y),fn=likelihood,C=C,invC=invC,detC=detC,x=x,method="L-BFGS-B",lower=c(10*.Machine$double.eps,rep(-Inf,tree$Nnode)),control=list(maxit=maxit))
	
	# return result
	states<-res$par[2:length(res$par)]
	names(states)<-c(length(tree$tip)+1,rownames(C)[(length(tree$tip)+1):nrow(C)])
	return(list(sig2=res$par[1],ace=states,logLik=-res$value,convergence=res$convergence))

}

