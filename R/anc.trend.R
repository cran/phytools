# this function estimates ancestral traits with a trend
# written by Liam J. Revell 2011

anc.trend<-function(phy,x,maxit=2000){

	# preliminaries
	# require dependencies
	require(ape)
	# set global
	tol<-1e-8
	# compute C
	D<-dist.nodes(phy)
	ntips<-length(phy$tip.label)
	Cii<-D[ntips+1,]
	C<-D; C[,]<-0
	counts<-vector()
	for(i in 1:nrow(D)) for(j in 1:ncol(D)) C[i,j]<-(Cii[i]+Cii[j]-D[i,j])/2
	dimnames(C)[[1]][1:length(phy$tip)]<-phy$tip.label
	dimnames(C)[[2]][1:length(phy$tip)]<-phy$tip.label
	C<-C[c(1:ntips,(ntips+2):nrow(C)),c(1:ntips,(ntips+2):ncol(C))]
	# sort x by phy$tip.label
	x<-x[phy$tip.label]

	# function returns the negative log-likelihood
	likelihood<-function(theta,x,C){
		a<-theta[1]
		u<-theta[2]
		sig2<-theta[3]
		y<-theta[4:length(theta)]
		logLik<-dmnorm(x=c(x,y),mean=(a+diag(C)*u),varcov=sig2*C,log=TRUE)
		return(-logLik)
	}

	# get reasonable starting values for the optimizer
	a<-mean(x)
	sig2<-var(x)/max(C)

	# perform ML optimization
	result<-optim(par=c(a,0,sig2,rep(a,phy$Nnode-1)),likelihood,x=x,C=C,method="L-BFGS-B",lower=c(-Inf,-Inf,tol,rep(-Inf,phy$Nnode-1)),control=list(maxit=maxit))

	# return the result
	ace<-c(result$par[c(1,4:length(result$par))]); names(ace)<-c(as.character(phy$edge[1,1]),rownames(C)[(length(phy$tip.label)+1):nrow(C)])
	return(list(ace=ace,mu=result$par[2],sig2=result$par[3],logL=-result$value,convergence=result$convergence,message=result$message))

}
