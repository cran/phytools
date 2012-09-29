# function computes phylogenetic variance-covariance matrix, including for internal nodes
# written by Liam J. Revell 2011

vcvPhylo<-function(tree,anc.nodes=T){
	n<-length(tree$tip.label)
	h<-nodeHeights(tree)[order(tree$edge[,2]),2]
	h<-c(h[1:n],0,h[(n+1):length(h)])
	M<-mrca(tree,full=anc.nodes)[c(1:n,anc.nodes*(n+2:tree$Nnode)),c(1:n,anc.nodes*(n+2:tree$Nnode))]
	C<-matrix(h[M],nrow(M),ncol(M))
	if(anc.nodes) rownames(C)<-colnames(C)<-c(tree$tip.label,n+2:tree$Nnode)
	else rownames(C)<-colnames(C)<-tree$tip.label
	return(C)
}

# returns the heights of each node
# written by Liam J. Revell 2011/2012

nodeHeights<-function(tree){
	if(attr(tree,"order")!="cladewise"||is.null(attr(tree,"order"))) t<-reorder(tree)
	else t<-tree
	root<-length(t$tip)+1
	X<-matrix(NA,nrow(t$edge),2)
	for(i in 1:nrow(t$edge)){
		if(t$edge[i,1]==root){
			X[i,1]<-0.0
			X[i,2]<-t$edge.length[i]
		} else {
			X[i,1]<-X[match(t$edge[i,1],t$edge[,2]),2]
			X[i,2]<-X[i,1]+t$edge.length[i]
		}
	}
	if(attr(tree,"order")!="cladewise"||is.null(attr(tree,"order")))
		o<-apply(matrix(tree$edge[,2]),1,function(x,y) which(x==y),y=t$edge[,2])
	else o<-1:nrow(t$edge)
	return(X[o,])
}
