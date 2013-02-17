# function adds a new tip to the tree
# written by Liam J. Revell 2012

bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL,position=0){
	if(is.null(where)) where<-length(tree$tip)+1
	if(is.null(edge.length)&&is.ultrametric(tree)){
		H<-nodeHeights(tree)
		if(where==(length(tree$tip)+1)) edge.length<-max(H)
		else edge.length<-max(H)-H[tree$edge[,2]==where,2]+position
	}
	tip<-list(edge=matrix(c(2,1),1,2),
		tip.label=tip.label,
		edge.length=edge.length,
		Nnode=1)
		class(tip)<-"phylo"
	obj<-bind.tree(tree,tip,where=where,position=position)
	return(obj)
}
