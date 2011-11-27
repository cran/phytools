# This function is a simplified REML version of brownie.lite()
# written by Liam J. Revell 2011

brownieREML<-function(tree,x,maxit=2000){
	# bookkeeping
	x<-as.matrix(x)
	n<-nrow(x) # number of species
	p<-ncol(tree$mapped.edge) # number of states

	# fit the single rate model
	likelihood1<-function(sig1,tree,x){
		phy<-scaleByMap(tree,rep(sig1,p))
		picX<-pic(x,phy,scaled=F,var.contrasts=T)
		logL<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
		return(-logL)
	}
	sig1<-mean(pic(x,tree)^2)
	logL1<--likelihood1(sig1,tree,x)

	# fit the multiple rate model
	likelihood2<-function(sig2,tree,x){
		phy<-scaleByMap(tree,sig2)
		picX<-pic(x,phy,scaled=F,var.contrasts=T)
		logL<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
		return(-logL)
	}
	res2<-optim(rep(1,p)*runif(n=p),likelihood2,tree=tree,x=x,method="L-BFGS-B",lower=rep(0,p))

	sig2<-res2$par; names(sig2)<-colnames(tree$mapped.edge)
	logL2<--res2$value
	convergence=(res2$convergence==0)

	return(list(sig2.single=sig1,logL1=logL1,sig2.multiple=sig2,logL2=logL2,convergence=convergence))

}

# This function scales a mapped tree by sig2
# written by Liam J. Revell 2011

scaleByMap<-function(mtree,sig2){
	edge.length<-rep(0,nrow(mtree$edge))
	for(i in 1:ncol(mtree$mapped.edge)) edge.length<-edge.length+sig2[i]*mtree$mapped.edge[,i]
	phy<-list(Nnode=mtree$Nnode,edge=mtree$edge,tip.label=mtree$tip.label,edge.length=edge.length)
	class(phy)<-"phylo"
	return(phy)
}
