# function to rescale simmap style trees
# written by Liam J. Revell 2012

rescaleSimmap<-function(tree,totalDepth=1.0){
	if(class(tree)=="multiPhylo"){
		trees<-list(); class(trees)<-"multiPhylo"
		trees<-lapply(tree,rescaleSimmap,totalDepth)
		return(trees)
	} else if(class(tree)=="phylo"){
		h<-max(nodeHeights(tree))
		s<-totalDepth/h
		tree$edge.length<-tree$edge.length*s
		maps<-lapply(tree$maps,"*",s)
		tree$maps<-maps
		tree$mapped.edge<-tree$mapped.edge*s
		return(tree)
	} else message("tree should be an object of class \"phylo\" or \"multiPhylo\"")
}
