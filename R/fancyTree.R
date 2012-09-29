# function to plot special types of phylogeny visualizations
# so far the implemented types are:
# "extinction" in which all branches leading to extinct taxa (or prior to the MRCA of extant species) are plotted with red dashed lines;
# "traitgram3d" which creates a 3D graph projecting the tree into two-dimensional morphospace (with time as the third axis)
# "droptip" creates a two panel plot with the tips to be pruned marked (panel 1) and then removed, and returns the pruned tree
# written by Liam J. Revell 2012

fancyTree<-function(tree,type=c("extinction","traitgram3d","droptip"),...,control=list()){
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	type<-matchType(type,c("extinction","traitgram3d","droptip"))
	if(type=="extinction") extinctionTree(tree)
	else if(type=="traitgram3d") traitgram3d(tree,...,control=control)
	else if(type=="droptip") return(droptipTree(tree,...))
	else stop(paste("do not recognize type = \"",type,"\"",sep=""))
}

# match type
# written by Liam J. Revell 2012

matchType<-function(type,types){
	for(i in 1:length(types))
		if(all(strsplit(type,split="")[[1]]==strsplit(types[i],split="")[[1]][1:length(strsplit(type,split="")[[1]])]))
			type=types[i]
	return(type)
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
		z<-setdiff(Descendants(tree,root.node,type="all"),Descendants(tree,ca,type="all"))
		edges[as.character(z)]<-1
	}
	z<-Descendants(tree,ca,type="all")
	y<-Descendants(tree,z)
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
		z<-setdiff(Descendants(tree,root.node,type="all"),Descendants(tree,ca,type="all"))
		edges[as.character(z)]<-1
	}
	z<-Descendants(tree,ca,type="all")
	y<-Descendants(tree,z)
	for(i in 1:length(z)) if(!any(tree$tip.label[y[[i]]]%in%keep)) edges[as.character(z[i])]<-1
	par(mfrow=c(2,1))
	plot.phylo(tree,edge.color=edges+1,edge.lty=edges+1,edge.width=2,no.margin=TRUE)
	dtree<-drop.tip(tree,tip); dtree$root.edge<-max(nodeHeights(tree))-max(nodeHeights(dtree))
	plot.phylo(dtree,edge.width=2,no.margin=TRUE,root.edge=TRUE)
	return(dtree)
}

