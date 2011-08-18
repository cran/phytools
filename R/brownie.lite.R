# This function fits two or more different evolutionary rates to different parts of a phylogenetic tree
# Takes a tree with a binary or multistate character "painted" on the branches in SIMMAP format 
# (see read.simmap()) and data for a single continuously valued trait
# Based on "Brownie" by O'Meara et al. 2006
# Written by Liam Revell 2010
brownie.lite<-function(tree,x,maxit=2000,test="chisq",nsim=100){
	# require dependencies
	require(ape)
	# bookkeeping
	x<-as.matrix(x)
	n<-nrow(x) # number of species
	p<-ncol(tree$mapped.edge) # number of states
	one<-matrix(1,n,1)
	# first compute C for the whole tree
	C1<-vcv.phylo(tree)
	C1<-C1[rownames(x),rownames(x)]
	# find the ancestral state vector
	a<-as.numeric(colSums(solve(C1))%*%x/sum(solve(C1)))
	# compute the MLE of sig1
	sig1<-as.numeric(t(x-one%*%a)%*%solve(C1)%*%(x-one%*%a)/n)
	# compute the log-likelihood
	likelihood.single<-function(theta,y,C){
		n<-length(y)
		one<-matrix(1,n,1)		
		sig<-theta[1]
		a<-theta[2]
		logL<-as.numeric(-t(x-one%*%t(a))%*%solve(sig*C)%*%(x-one%*%t(a))/2-n*log(2*pi)/2-determinant(sig*C)$modulus[1]/2)
		return(-logL)
	}
	# optimize single rate model (we can feed it our ML estimates here)
	model1<-optim(c(sig1,a),fn=likelihood.single,y=x,C=C1,control=list(maxit=maxit),hessian=TRUE,method="L-BFGS-B",lower=c(0,-Inf))
	# convert for convenience of output
	sig1<-model1$par[1]
	logL1<--model1$value
	vcv1<-solve(model1$hessian); rownames(vcv1)<-c("sig","a"); colnames(vcv1)<-rownames(vcv1)
	a1<-model1$par[2]
	# compute separate C for each state
	multi.tre<-list(); class(multi.tre)<-"multiPhylo"
	C2<-array(dim=c(nrow(C1),ncol(C1),ncol(tree$mapped.edge)))
	for(i in 1:ncol(tree$mapped.edge)){
		multi.tre[[i]]<-tree
		multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		C2[,,i]<-temp[rownames(x),rownames(x)]
	}
	# compute the log-likelihood
	likelihood.multiple<-function(theta,y,C){
		n<-length(y); p<-dim(C)[3]
		one<-matrix(1,n,1)
		sig<-vector(mode="numeric")
		for(i in 1:p) sig[i]<-theta[i]
		a<-theta[p+1]
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+sig[i]*C[,,i]
		logL<--t(y-one%*%a)%*%solve(V)%*%(y-one%*%a)/2-n*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}
	# compute the starting parameter values
	starting<-vector(mode="numeric")
	lower<-vector(mode="numeric")
	for(i in 1:p){ 
		starting[i]<-sig1
		lower[i]<-0.001*sig1 # this sets the lower limit on the sigmas during optimization to be 0.1% of sigma in the single rate model
	}
	starting[p+1]<-a; lower[p+1]<--Inf
	# optimize using built-in optimizer
	model2<-optim(starting,fn=likelihood.multiple,y=x,C=C2,control=list(maxit=maxit),hessian=TRUE,method="L-BFGS-B",lower=lower)
	# log-likelihood for the multi-matrix model
	logL2<--model2$value
	while(logL2<logL1){
		message("False convergence on first try; trying again with new starting values.")
		model2<-optim(starting*2*runif(n=length(starting)),fn=likelihood.multiple,y=x,C=C2,control=list(maxit=maxit),hessian=TRUE,method="L-BFGS-B",lower=lower)
		logL2<--model2$value
	}
	# convert for convenience of output
	sig.i<-vector(mode="numeric"); states<-vector(mode="character")
	for(i in 1:p){
		sig.i[i]<-model2$par[i]
		states[i]<-multi.tre[[i]]$state
	}
	names(sig.i)<-states
	a2<-model2$par[p+1]
	# covariance matrix
	vcv2<-solve(model2$hessian); rownames(vcv2)<-c(states,"a"); colnames(vcv2)<-rownames(vcv2)
	# report convergence
	if(model2$convergence==0){
		converged<-"Optimization has converged."
	} else {
		converged<-"Optimization may not have converged.  Consider increasing maxit."
	}
	# return results
	if(test=="chisq"){
		return(list(sig2.single=sig1,a.single=a1,var.single=vcv1[1,1],logL1=logL1,k1=2,sig2.multiple=sig.i,a.multiple=a2,vcv.multiple=vcv2[1:p,1:p],logL.multiple=logL2,k2=p+1,P.chisq=pchisq(2*(logL2-as.numeric(logL1)),p-1,lower.tail=FALSE),convergence=converged))
	} else if(test=="simulation"){
		# compute the observed LR
		LR<-2*(logL2-as.numeric(logL1))
		# simulate some data under the null
		X<-fastBM(tree,a=a1,sig2=sig1,nsim=nsim)
		# now compute the P-value based on simulation
		Psim<-0
		for(i in 1:nsim){
			sim<-brownie.lite(tree,X[,i])
			while(sim$convergence!="Optimization has converged."||(sim$logL.multiple-sim$logL1)<0){
				X[,i]<-fastBM(tree,a=a1,sig2=sig1)
				sim<-brownie.lite(tree,X[,i])
			}
			Psim<-Psim+(LR>2*(sim$logL.multiple-sim$logL1))/nsim
		}
		return(list(sig2.single=sig1,a.single=a1,var.single=vcv1[1,1],logL1=logL1,k1=2,sig2.multiple=sig.i,a.multiple=a2,vcv.multiple=vcv2[1:p,1:p],logL.multiple=logL2,k2=p+1,P.sim=Psim,convergence=converged))
	} else {
		print(paste("test =",test,"is not a valid option - usingn chi-sq test"))
		return(list(sig2.single=sig1,a.single=a1,var.single=vcv1[1,1],logL1=logL1,k1=2,sig2.multiple=sig.i,a.multiple=a2,vcv.multiple=vcv2[1:p,1:p],logL.multiple=logL2,k2=p+1,P.chisq=pchisq(2*(logL2-as.numeric(logL1)),p-1,lower.tail=FALSE),convergence=converged))
	}		
}
