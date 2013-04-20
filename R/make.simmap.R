# function creates a stochastic character mapped tree as a modified "phylo" object
# written by Liam Revell 2013

make.simmap<-function(tree,x,model="SYM",nsim=1,...){
	if(class(tree)=="multiPhylo"){
		ff<-function(yy,x,model,nsim,...){
			zz<-make.simmap(yy,x,model,nsim,...)
			if(nsim>1) class(zz)<-NULL
			return(zz)
		}	
		if(nsim>1) mtrees<-unlist(sapply(tree,ff,x,model,nsim,...,simplify=FALSE),recursive=FALSE)
		else mtrees<-sapply(tree,ff,x,model,nsim,...,simplify=FALSE)
		class(mtrees)<-"multiPhylo"
	} else {
		if(hasArg(pi)) pi<-list(...)$pi
		else pi<-"equal"
		if(hasArg(message)) pm<-list(...)$message
		else pm<-TRUE
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-1e-8
		# check
		if(class(tree)!="phylo") stop("'tree' should be an object of class 'phylo'")
		# if vector convert to binary matrix
		if(!is.matrix(x)) xx<-to.matrix(x,sort(unique(x)))
		else xx<-x
		xx<-xx[tree$tip.label,]
		xx<-xx/rowSums(xx)
		# reorder to cladewise
		tree<-bt<-reorder.phylo(tree,"cladewise")
		if(!is.binary.tree(bt)) bt<-multi2di(bt)
		# some preliminaries
		N<-length(tree$tip)
		m<-ncol(xx)
		root<-N+1
		# get conditional likelihoods & model
		XX<-apeAce(bt,xx,model)
		II<-XX$index.matrix
		L<-XX$lik.anc; rownames(L)<-length(bt$tip)+1:nrow(L)
		if(!is.binary.tree(tree)){
			ancNames<-matchNodes(tree,bt)
			L<-L[as.character(ancNames[,2]),]
			rownames(L)<-ancNames[,1]
		}
		L<-rbind(xx,L)
		rownames(L)[1:N]<-1:N
		if(any(XX$rates<tol)){
			message(paste("\nWarning: some elements of Q not numerically distinct from 0; setting to",tol,"\n"))
			XX$rates[XX$rates<tol]<-tol
		}
		Q<-matrix(XX$rates[II],m,m,dimnames=list(colnames(L),colnames(L)))
		diag(Q)<--rowSums(Q,na.rm=TRUE)
		if(pi[1]=="equal") pi<-setNames(rep(1/m,m),colnames(L)) # set equal
		else if(pi[1]=="estimated") pi<-statdist(Q) # set from stationary distribution
		else pi<-pi/sum(pi) # obtain from input
		if(pm) printmessage(Q,pi)
		mtrees<-replicate(nsim,smap(tree,x,N,m,root,L,Q,pi),simplify=FALSE)
		if(length(mtrees)==1) mtrees<-mtrees[[1]]
		else class(mtrees)<-"multiPhylo"
	}	
	(if(hasArg(message)) list(...)$message else TRUE)
	if((if(hasArg(message)) list(...)$message else TRUE)&&class(tree)=="phylo") message("Done.")
	return(mtrees)
}

# convert vector of x to binary matrix
# written by Liam J. Revell 2012
to.matrix<-function(x,seq){
	X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
	for(i in 1:length(seq)) X[x==seq[i],i]<-1
	return(X)
}

# function does the stochastic mapping, conditioned on our model & given the conditional likelihoods
# written by Liam J. Revell 2013
smap<-function(tree,x,N,m,root,L,Q,pi){
	# create the map tree object
	mtree<-tree; mtree$maps<-list()
	mtree$mapped.edge<-matrix(0,nrow(tree$edge),m,dimnames=list(paste(tree$edge[,1],",",tree$edge[,2],sep=""),colnames(L)))
	# now we want to simulate the node states & histories by pre-order traversal
	NN<-matrix(NA,nrow(tree$edge),2) # our node values
	NN[which(tree$edge[,1]==root),1]<-rstate(L[as.character(root),]*pi/sum(L[as.character(root),]*pi)) # assign root
	for(j in 1:nrow(tree$edge)){
		# conditioned on the start value, assign end value of node (if internal)
		p<-expm(Q*tree$edge.length[j])[NN[j,1],]*L[as.character(tree$edge[j,2]),]
		NN[j,2]<-rstate(p/sum(p))
		NN[which(tree$edge[,1]==tree$edge[j,2]),1]<-NN[j,2]
		# now simulate on the branches
		accept<-FALSE	
		while(!accept){
			map<-sch(NN[j,1],tree$edge.length[j],Q)
			if(names(map)[length(map)]==NN[j,2]) accept=TRUE
		}	
		mtree$maps[[j]]<-map
		for(k in 1:length(mtree$maps[[j]])) 
			mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
	}
	return(mtree)
}

# function generates a history along a branch
# written by Liam J. Revell 2013
sch<-function(start,t,Q){
	dt<-setNames(0,start)
	while(sum(dt)<t){
		s<-names(dt)[length(dt)]
		dt[length(dt)]<-rexp(n=1,rate=-Q[s,s])
		if(sum(dt)<t){
			dt<-c(dt,0)
			names(dt)[length(dt)]<-rstate(Q[,s][-match(s,rownames(Q))]/sum(Q[,s][-match(s,rownames(Q))]))
		} else dt[length(dt)]<-dt[length(dt)]-sum(dt)+t
	}
	return(dt)
}

# function returns random state with probability given by y
# written by Liam J. Revell 2013
rstate<-function(y){
	if(length(y)==1) return(names(y)[1])
	else return(names(which(rmultinom(1,1,y/sum(y))[,1]==1)))
}

# function uses numerical optimization to solve for the stationary distribution
# written by Liam J. Revell 2013
statdist<-function(Q){
	foo<-function(theta,Q){
		Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
		sum((Pi%*%Q)^2)
	}
	k<-nrow(Q)
	if(nrow(Q)>2){ 
		fit<-optim(rep(1/k,k-1),foo,Q=Q,control=list(reltol=1e-16))
		return(setNames(c(fit$par[1:(k-1)],1-sum(fit$par[1:(k-1)])),rownames(Q)))
	} else {
		fit<-optimize(foo,interval=c(0,1),Q=Q)
		return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
	}
}

# print message
# written by Liam J. Revell 2013
printmessage<-function(Q,pi){
	message("make.simmap is sampling character histories conditioned on the transition matrix\nQ =")
	print(Q)
	message("(estimated using likelihood);")
	message("and root node prior probabilities\npi =")
	print(pi)
}

# function for conditional likelihoods at nodes, from ace(...,type="discrete")
# modified (only very slightly) from E. Paradis et al. 2013
apeAce<-function(tree,x,model){
	ip<-0.1
	nb.tip<-length(tree$tip.label)
	nb.node<-tree$Nnode
	if(is.matrix(x)){ 
		x<-x[tree$tip.label,]
		nl<-ncol(x)
		lvls<-colnames(x)
	} else {
		x<-x[tree$tip.label]
  		if(!is.factor(x)) x<-factor(x)
		nl<-nlevels(x)
		lvls<-levels(x)
		x<-as.integer(x)
	}
	if(is.character(model)){
		rate<-matrix(NA,nl,nl)
		if(model=="ER") np<-rate[]<-1
		if(model=="ARD"){
			np<-nl*(nl-1)
			rate[col(rate)!=row(rate)]<-1:np
		}
		if (model=="SYM") {
			np<-nl*(nl-1)/2
			sel<-col(rate)<row(rate)
			rate[sel]<-1:np
			rate<-t(rate)
			rate[sel]<-1:np
		}
	} else {
		if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
		if(ncol(model)!=nl) stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
		rate<-model
		np<-max(rate)
	}
	index.matrix<-rate
	tmp<-cbind(1:nl,1:nl)
	index.matrix[tmp]<-NA
	rate[tmp]<-0
	rate[rate==0]<-np+1
	liks<-matrix(0,nb.tip+nb.node,nl)
	TIPS<-1:nb.tip
	if(is.matrix(x)) liks[TIPS,]<-x
	else liks[cbind(TIPS,x)]<-1
	phy<-reorder(tree,"pruningwise")
	Q<-matrix(0,nl,nl)
	dev<-function(p,output.liks=FALSE){
		if(any(is.nan(p))||any(is.infinite(p))) return(1e50)
		comp<-numeric(nb.tip+nb.node)
		Q[]<-c(p,0)[rate]
		diag(Q)<--rowSums(Q)
		for(i in seq(from=1,by=2,length.out=nb.node)){
			j<-i+1L
			anc<-phy$edge[i,1]
			des1<-phy$edge[i,2]
			des2<-phy$edge[j,2]
			v.l<-matexpo(Q*phy$edge.length[i])%*%liks[des1,]
			v.r<-matexpo(Q*phy$edge.length[j])%*%liks[des2,]
			v<-v.l*v.r
			comp[anc]<-sum(v)
			liks[anc,]<-v/comp[anc]
		}
		if(output.liks) return(liks[-TIPS,])
		dev<--2*sum(log(comp[-TIPS]))
		if(is.na(dev)) Inf else dev
	}
	out<-nlminb(rep(ip,length.out=np),function(p) dev(p),lower=rep(0,np),upper=rep(1e50,np))
	obj<-list()
	obj$loglik<--out$objective/2
	obj$rates<-out$par
	obj$index.matrix<-index.matrix
	obj$lik.anc<-dev(obj$rates,TRUE)
	colnames(obj$lik.anc)<-lvls
	return(obj)
}
