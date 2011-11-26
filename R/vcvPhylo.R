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
