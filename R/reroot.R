# function to re-root a phylogeny along an edge
# written by Liam Revell 2011, 2013

reroot<-function(tree,node.number,position){
	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	if(node.number>length(tree$tip)){
		# first, re-root the tree at node.number
		tr<-root(tree,node=node.number,resolve.root=T)
		# now re-allocate branch length to the two edges descending from the new root node
		b<-sum(tr$edge.length[tr$edge==(length(tree$tip)+1)])
		tr$edge.length[tr$edge==(length(tree$tip)+1)]<-c(position,b-position)
	} else {
		# first, root the tree at the parent of node.number
		tr1<-root(tree,node=tree$edge[match(node.number,tree$edge[,2]),1])
		# now drop tip
		tr1<-drop.tip(tr1,tree$tip.label[node.number])
		# create phylo object
		tr2<-list(edge=matrix(c(3L,1L,3L,2L),2,2,byrow=T),tip.label=c(tree$tip.label[node.number],"NA"),edge.length=c(tree$edge.length[match(node.number,tree$edge[,2])]-position,position),Nnode=1)
		class(tr2)<-"phylo"
		tr<-bind.tree(tr2,tr1,where=which(tr2$tip.label=="NA"))
	}
	return(tr)
}
