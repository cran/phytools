# function plots a tree
# written by Liam Revell 2011

plotTree<-function(tree,fsize=1.0,ftype="reg",lwd=2,pts=TRUE){
	# check font
	ftype<-which(c("reg","b","i","bi")==ftype)
	# check tree
	if(class(tree)!="phylo") stop("tree should be object of class 'phylo.'")
	# reorder
	cw<-reorder(tree)
	pw<-reorder(tree,"pruningwise")
	# count nodes and tips
	n<-length(cw$tip); m<-cw$Nnode
	# Y coordinates for nodes
	Y<-matrix(NA,m+n,1)
	# first, assign y coordinates to all the tip nodes
	Y[cw$edge[cw$edge[,2]<=length(cw$tip),2]]<-1:n
	# get Y coordinates of the nodes
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
		Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
	}
	# compute node heights
	root<-length(cw$tip)+1
	node.height<-matrix(NA,nrow(cw$edge),2)
	for(i in 1:nrow(cw$edge)){
		if(cw$edge[i,1]==root){
			node.height[i,1]<-0.0
			node.height[i,2]<-cw$edge.length[i]
		} else {
			node.height[i,1]<-node.height[match(cw$edge[i,1],cw$edge[,2]),2]
			node.height[i,2]<-node.height[i,1]+cw$edge.length[i]
		}
	}
	# open plot
	par(mar=c(0.1,0.1,0.1,0.1))
	plot.new(); plot.window(xlim=c(0,1.05*(max(node.height)+fsize*max(strwidth(cw$tip.label)))),ylim=c(1,max(Y)))
	for(i in 1:nrow(cw$edge)){
		if(pts) points(node.height[i,],c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),pch=20)
 		lines(node.height[i,],c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),lwd=lwd,lend=0)
	}
	for(i in 1:m) lines(node.height[which(cw$edge[,1]==nodes[i]),1],Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],lwd=lwd)
	for(i in 1:n) text(node.height[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=4,cex=fsize,font=ftype)
}
