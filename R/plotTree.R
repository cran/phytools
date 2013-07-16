# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012, 2013

plotTree<-function(tree,...){
	if(hasArg(color)) color<-list(...)$color
	else color<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-"reg"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(pts)) pts<-list(...)$pts
	else pts<-FALSE
	if(hasArg(node.numbers)) node.numbers<-list(...)$node.numbers
	else node.numbers<-FALSE
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-NULL
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(setEnv)) setEnv<-list(...)$setEnv
	else setEnv<-FALSE
	if(hasArg(part)) part<-list(...)$part
	else part<-1.0
	if(class(tree)=="multiPhylo"){
		par(ask=TRUE)
		if(!is.null(color)) names(color)<-"1"
		for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,type=type,setEnv=setEnv,part=part)
	} else {
		tree$maps<-as.list(tree$edge.length)
		for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
		if(!is.null(color)) names(color)<-"1"
		plotSimmap(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,type=type,setEnv=setEnv,part=part)
	}
}
