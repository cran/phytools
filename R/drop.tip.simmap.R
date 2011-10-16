# this function drops tips from the modified "phylo" format trees created by read.simmap() and maintains the $mapped.edge matrix
# written by Liam Revell, 2010
drop.tip.simmap<-function(tree,tip){
	# set tolerance
	tol<-1e-8
	# create the basic tree structure
	new.tree<-drop.tip(tree,tip)
	new.tree$mapped.edge<-matrix(NA,nrow(new.tree$edge),ncol(tree$mapped.edge),dimnames=list(edge=apply(new.tree$edge,1,function(x) paste(x,collapse=",")),state=colnames(tree$mapped.edge)))
	j<-1
	for(i in 1:ncol(tree$mapped.edge)){
		phy<-tree
		phy$edge.length<-phy$mapped.edge[,i]
		phy<-drop.tip(phy,tip)
		new.tree$mapped.edge[,j]<-phy$edge.length
		colnames(new.tree$mapped.edge)[j]<-colnames(phy$mapped.edge)[i]
		if(sum(new.tree$mapped.edge[,j])>tol) j<-j+1 # this has the effect of eliminating any mapped traits absent from the pruned tree
	}
	new.tree$mapped.edge<-matrix(new.tree$mapped.edge[,1:(j-1)],nrow(new.tree$mapped.edge),j-1,dimnames=list(edge=rownames(new.tree$mapped.edge),state=colnames(new.tree$mapped.edge)[1:(j-1)]))
	return(new.tree)
}
