# function slices a tree at slice and returns all subtrees
# it uses extract.clade(), if trivial==FALSE subtrees with length than 2 taxa are ignored
# written by Liam Revell 2011/2012

treeSlice<-function(tree,slice,trivial=FALSE){
	if(class(tree)!="phylo") stop("tree should be object of class 'phylo.'")
	tree<-reorder(tree) # reorder cladewise
	H<-nodeHeights(tree)
	edges<-which(H[,2]>slice&H[,1]<slice)
	nodes<-tree$edge[edges,2]
	if(!trivial) nodes<-nodes[nodes>length(tree$tip)]
	trees<-list(); class(trees)<-"multiPhylo"
	for(i in 1:length(nodes)){
		if(nodes[i]>length(tree$tip)){ 
			trees[[i]]<-extract.clade(tree,node=nodes[i])
			trees[[i]]$root.edge<-H[which(tree$edge[,2]==nodes[i]),2]-slice
		} else {
			z<-list(edge=matrix(c(2,1),1,2),edge.length=H[which(tree$edge[,2]==nodes[i]),2]-slice,tip.label=tree$tip.label[nodes[i]],Nnode=1L)
			class(z)<-"phylo"; trees[[i]]<-z
		}
	}
	return(trees)
}

