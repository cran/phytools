# this function fits a "diversity-dependent-evolutionary-diversification" model (similar to Mahler et al. 2010)
# written by Liam Revell, 2010/2011

fitDiversityModel<-function(tree,x,d=NULL,showTree=TRUE){
	# check for & load "ape"	
	if(!require(ape)) stop("must first install 'ape' package.") # require ape	
	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	if(!is.binary.tree(tree)){
		message("tree must be fully bifurcating. randomly resolving.")
		tree<-multi2di(tree)
	}
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)) x<-x[,1]
	if(is.null(names(x))){
		if(length(x)==length(tree$tip)){
			message("x has no names; assuming x is in the same order as tree$tip.label")
			names(x)<-tree$tip.label
		} else
			stop("x has no names and is a different length than tree$tip.label")
	}
	if(any(is.na(match(tree$tip.label,names(x))))){
		message("some species in tree are missing from data, dropping missing taxa from the tree")
		tree<-drop.tip(tree,tree$tip.label[-match(names(x),tree$tip.label)])
	}
	if(any(is.na(match(names(x),tree$tip.label)))){
		message("some species in data are missing from tree, dropping missing taxa from the data")
		x<-x[tree$tip.label]
	}
	if(any(is.na(x))){
		message("some data given as 'NA', dropping corresponding species from tree")
		tree<-drop.tip(tree,names(which(is.na(x))))
	}
	if(!is.null(d)){
		if(is.data.frame(d)) d<-as.matrix(d)
		if(is.matrix(d)) d<-d[,1]
		if(is.null(names(d))){
			if(length(d)==tree$Nnode){
				message("d has no names; assuming d is in node number order of the resolved tree")
				names(d)<-c(length(tree$tip)+1:tree$Nnode)
			} else
				stop("d has no names and is a different length than tree$Nnode for the resolved tree")
		}
	}
	# tolerance
	tol<-1e-8
	# compute tree length
	tree.length<-mean(diag(vcv(tree)))
	# compute contrasts
	pic.x<-pic(x,tree)
	# if computing d
	if(is.null(d)){
		# compute lineage diversity at each node
		ages<-branching.times(tree)
		d<-vector()
		for(i in 1:length(ages)) d[i]<-sum(ages>ages[i])
		names(d)<-names(ages)
	}
	maxd<-max(d)
	d<-d/(maxd+1)
	# likelihood function
	likelihood<-function(theta,y,phy,diversity){
		sig0<-theta[1]
		scaled.psi<-theta[2]
		for(i in 1:nrow(phy$edge)){
			vi<-phy$edge.length[i]
			phy$edge.length[i]<-vi+vi*scaled.psi*diversity[as.character(phy$edge[i,1])]
		}
		D<-vcv(phy)*sig0
		D<-D[names(y),names(y)]
		Dinv<-solve(D)
		a<-as.numeric(colSums(Dinv)%*%y/sum(Dinv))
		logL<-as.numeric(-t(y-a)%*%Dinv%*%(y-a)/2-determinant(D)$modulus[1]/2-length(y)*log(2*pi)/2)
		if(showTree) plot(phy)
		return(-logL)
	}
	# optimize
	res=optim(c(mean(pic.x^2),0),likelihood,y=x,phy=tree,diversity=d,method="L-BFGS-B",lower=c(tol,-1),upper=c(Inf,1),hessian=TRUE)	
	# return result
	if(var(d)>0)
		return(list(logL=-res$value,sig0=res$par[1],psi=res$par[2]*res$par[1]/(maxd+1),vcv=matrix(solve(res$hessian),2,2,dimnames=list(c("sig0","psi"),c("sig0","psi"))),convergence=(res$convergence==0)))
	else {
		message("psi not estimable because diversity is constant through time.")
		return(list(logL=-res$value,sig0=res$par[1],vcv=matrix(solve(res$hessian[1,1]),1,1,dimnames=list(c("sig0"),c("sig0"))),convergence=(res$convergence==0)))
	}
}
