# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012

plotTree<-function(tree,color=NULL,fsize=1.0,ftype="reg",lwd=2,pts=TRUE,node.numbers=FALSE){
	if(class(tree)=="multiPhylo"){
		par(ask=TRUE)
		if(!is.null(color)) names(color)<-"1"
		for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers)
	} else {
		tree$maps<-as.list(tree$edge.length)
		for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
		if(!is.null(color)) names(color)<-"1"
		plotSimmap(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers)
	}
}
