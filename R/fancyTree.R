# function to plot special types of phylogeny visualizations
# so far the implemented types are:
# "extinction" in which all branches leading to extinct taxa (or prior to the MRCA of extant species) are plotted with red dashed lines;
# "traitgram3d" which creates a 3D graph projecting the tree into two-dimensional morphospace (with time as the third axis)
# "droptip" creates a two panel plot with the tips to be pruned marked (panel 1) and then removed, and returns the pruned tree
# "xkcd" creates an xkcd-comic style phylogeny
# "densitymap" maps the posterior density of a binary stochastic character mapping
# "contmap" maps reconstructed trait evolution for a continuous character on the tree
# written by Liam J. Revell 2012, 2013

fancyTree<-function(tree,type=c("extinction","traitgram3d","droptip","xkcd","densitymap","contmap"),...,control=list()){
	type<-matchType(type,c("extinction","traitgram3d","droptip","xkcd","densitymap"))
	if(class(tree)!="phylo"&&type%in%c("extinction","traitgram3d","droptip","xkcd")) stop("tree should be an object of class 'phylo'")
	else if(class(tree)!="multiPhylo"&&type=="densitymap") stop("for type='densitymap' tree should be an object of class 'multiPhylo'")
	if(type=="extinction") extinctionTree(tree)
	else if(type=="traitgram3d") traitgram3d(tree,...,control=control)
	else if(type=="droptip") return(droptipTree(tree,...))
	else if(type=="xkcd") plotXkcdTree(tree,...)
	else if(type=="densitymap") plotDensityMap(tree,...)
	else if(type=="contmap") plotContMap(tree,...)
	else stop(paste("do not recognize type = \"",type,"\"",sep=""))
}

# extinctionTree internal function
# written by Liam J. Revell 2012

extinctionTree<-function(tree){
	edges<-rep(0,nrow(tree$edge))
	names(edges)<-tree$edge[,2]
	extant<-getExtant(tree)
	ca<-findMRCA(tree,extant)
	root.node<-length(tree$tip)+1
	if(ca!=root.node){
		z<-setdiff(getDescendants(tree,root.node),getDescendants(tree,ca))
		edges[as.character(z)]<-1
	}
	z<-getDescendants(tree,ca)
	y<-lapply(z,getDescendants,tree=tree)
	for(i in 1:length(z)) if(!any(tree$tip.label[y[[i]]]%in%extant)) edges[as.character(z[i])]<-1
	plot.phylo(tree,edge.color=edges+1,edge.lty=edges+1,edge.width=2,no.margin=TRUE)
}

# traitgram3d internal function
# written by Liam J. Revell 2012

traitgram3d<-function(tree,...,control){
	if(hasArg(X)) X<-list(...)$X
	else stop("no phenotypic data provided")
	if(!hasArg(A)){
		if(is.null(control$maxit)) maxit<-2000
		else maxit<-control$maxit
		Y<-apply(X,2,function(x,tree) anc.ML(tree,x,maxit),tree=tree)
		convergence<-sapply(Y,function(x) x$convergence)
		if(any(convergence!=0)) warning("anc.ML may not have converged; consider increasing maxit.")
		A<-sapply(Y,function(x) x$ace)
	} else { 
		A<-list(...)$A
		A<-A[as.character(1:tree$Nnode+length(tree$tip)),]
	}
	if(is.null(colnames(X))) colnames(X)<-c("x","y")
	X<-cbind(X,diag(vcv(tree))[rownames(X)])
	A<-cbind(A,nodeHeights(tree)[match(rownames(A)[1:nrow(A)],tree$edge)])
	colnames(X)[3]<-colnames(A)[3]<-"time"
	phylomorphospace3d(tree,X,A,control=control)
}

# droptipTree internal function
# written by Liam J. Revell 2012

droptipTree<-function(tree,...){
	if(hasArg(tip)) tip<-list(...)$tip
	else stop("need to provide tip or tips to drop")
	edges<-rep(0,nrow(tree$edge))
	names(edges)<-tree$edge[,2]
	keep<-setdiff(tree$tip.label,tip)
	ca<-findMRCA(tree,keep)
	root.node<-length(tree$tip)+1
	if(ca!=root.node){
		z<-setdiff(getDescendants(tree,root.node),getDescendants(tree,ca))
		edges[as.character(z)]<-1
	}
	z<-getDescendants(tree,ca)
	foo<-function(x,tree){
		n<-length(tree$tip.label)
		y<-getDescendants(tree,x)
		y<-y[y<=n]
		return(y)
	}
	y<-lapply(z,foo,tree=tree)
	for(i in 1:length(z)) if(!any(tree$tip.label[y[[i]]]%in%keep)) edges[as.character(z[i])]<-1
	par(mfrow=c(2,1))
	plot.phylo(tree,edge.color=edges+1,edge.lty=edges+1,edge.width=2,no.margin=TRUE)
	dtree<-drop.tip(tree,tip); dtree$root.edge<-max(nodeHeights(tree))-max(nodeHeights(dtree))
	plot.phylo(dtree,edge.width=2,no.margin=TRUE,root.edge=TRUE)
	return(dtree)
}

# plotXkcdTree internal function
# written by Liam J. Revell 2012

plotXkcdTree<-function(tree,...){
	if(hasArg(file)) file<-list(...)$file
	else file<-NULL
	if(hasArg(gsPath)) gsPath<-list(...)$gsPath
	else gsPath<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-2
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-7
	if(hasArg(color)) color<-list(...)$color
	else color<-"blue"
	if(hasArg(dim)) dim<-list(...)$dim
	else dim<-c(8.5,11)
	if(hasArg(jitter)) jitter<-list(...)$jitter
	else jitter<-0.001
	if(hasArg(waver)) waver<-list(...)$waver
	else waver<-c(0.1,0.1)
	if(hasArg(tilt)) tilt<-list(...)$tilt
	else tilt<-0
	if(hasArg(right)) right<-list(...)$right
	else right<-TRUE
	xkcdTree(tree,file,gsPath,fsize,lwd,color,dim,jitter,waver,tilt,right)
}

# plotDensityMap internal function
# written by Liam J. Revell 2012

plotDensityMap<-function(trees,...){
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-3
	if(hasArg(check)) check<-list(...)$check
	else check<-FALSE
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-FALSE
	densityMap(trees,res,fsize,check,legend,outline)
}

# plotContMap internal function
# written by Liam J. Revell 2012

plotContMap<-function(tree,...){
	if(hasArg(x)) x<-list(...)$x
	else stop("need to provide vector 'x' of phenotypic trait values")
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-4
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(lims)) lims<-list(...)$lims
	else lims<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-TRUE
	if(hasArg(sig)) sig<-list(...)$sig
	else sig<-3
	contMap(tree,x,res,fsize,ftype,lwd,legend,lims,outline,sig)
}


