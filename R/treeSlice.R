# function slices a tree at slice and returns all subtrees
# it uses extract.clade() and ignores subtrees with fewer than 2 tips
# written by Liam Revell 2011

treeSlice<-function(tree,slice){
	if(class(tree)!="phylo") stop("tree should be object of class 'phylo.'")
	tree<-reorder(tree) # reorder cladewise
	# compute node heights
	root<-length(tree$tip)+1
	node.height<-matrix(NA,nrow(tree$edge),2)
	for(i in 1:nrow(tree$edge)){
		if(tree$edge[i,1]==root){
			node.height[i,1]<-0.0
			node.height[i,2]<-tree$edge.length[i]
		} else {
			node.height[i,1]<-node.height[match(tree$edge[i,1],tree$edge[,2]),2]
			node.height[i,2]<-node.height[i,1]+tree$edge.length[i]
		}
	}
	edges<-which(node.height[,2]>slice&node.height[,1]<slice)
	nodes<-tree$edge[edges,2]; nodes<-nodes[nodes>length(tree$tip)]
	trees<-list(); class(trees)<-"multiPhylo"
	for(i in 1:length(nodes)){ 
		trees[[i]]<-extract.clade(tree,node=nodes[i])
		trees[[i]]$root.edge<-node.height[which(tree$edge[,2]==nodes[i]),2]-slice
	}
	return(trees)
}

