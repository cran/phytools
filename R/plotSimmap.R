# function plots a stochastic character mapped tree
# written by Liam Revell 2011

plotSimmap<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,pts=TRUE,node.numbers=FALSE){
	if(class(tree)=="multiPhylo"){
		par(ask=TRUE)
		for(i in 1:length(tree)) plotSimmap(tree[[i]],colors=colors,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers)
	} else {
		# check font
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!ftype) fsize=0 
		# check colors
		if(is.null(colors)){ 
			colors<-palette()
			names(colors)<-as.character(1:8)
		}
		# check tree
		if(class(tree)!="phylo") stop("tree should be object of class 'phylo.'")
		if(is.null(tree$maps)) stop("tree should contain mapped states on edges.")
		# swap out "_" character for spaces (assumes _ is a place holder)
		tree$tip.label<-gsub("_"," ",tree$tip.label)
		# reorder
		cw<-reorderSimmap(tree)
		pw<-reorderSimmap(tree,"pruningwise")
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
		plot.new()
		if(fsize*max(strwidth(cw$tip.label))<1.0){
			c<-(1-fsize*max(strwidth(cw$tip.label)))/max(node.height)
			cw$edge.length<-c*cw$edge.length
			cw$maps<-lapply(cw$maps,function(x) x<-c*x)
			node.height<-c*node.height
		} else message("Font size too large to properly rescale tree to window.")
		plot.window(xlim=c(0,max(node.height)+fsize*max(strwidth(cw$tip.label))),ylim=c(1,max(Y)))
		for(i in 1:m) lines(node.height[which(cw$edge[,1]==nodes[i]),1],Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],col=colors[names(cw$maps[[match(nodes[i],cw$edge[,1])]])[1]],lwd=lwd)
		for(i in 1:nrow(cw$edge)){
			x<-node.height[i,1]
	 		for(j in 1:length(cw$maps[[i]])){
				lines(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=2)
				if(pts) points(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),pch=20,lwd=(lwd-1))
				x<-x+cw$maps[[i]][j]; j<-j+1
			}
		}
		if(node.numbers){
			symbols(0,mean(Y[cw$edge[cw$edge[,1]==(length(cw$tip)+1),2]]),rectangles=matrix(c(1.2*fsize*strwidth(as.character(length(cw$tip)+1)),1.4*fsize*strheight(as.character(length(cw$tip)+1))),1,2),inches=F,bg="white",add=T)
			text(0,mean(Y[cw$edge[cw$edge[,1]==(length(cw$tip)+1),2]]),length(cw$tip)+1,cex=fsize)
			for(i in 1:nrow(cw$edge)){
				x<-node.height[i,2]
				if(cw$edge[i,2]>length(tree$tip)){
					symbols(x,Y[cw$edge[i,2]],rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=F,bg="white",add=T)
					text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
				}
			}
		}
		for(i in 1:n) if(ftype) text(node.height[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=4,offset=0.2*lwd/3+0.2/3,cex=fsize,font=ftype)
	}
	# reset margin
	par(mar=c(5,4,4,2)+0.1)
}
	
# function reorders simmap tree
# written Liam Revell 2011
	
reorderSimmap<-function(tree,order="cladewise"){
	ntree<-reorder(tree,order)
	o<-whichorder(ntree$edge[,2],tree$edge[,2])
	ntree$mapped.edge<-tree$mapped.edge[o,]
	ntree$maps<-tree$maps[o]
	return(ntree)
}

# function whichorder
# written by Liam Revell 2011

whichorder<-function(x,y){
	n<-length(x); ind<-vector()
	for(i in 1:n) ind[i]<-which(x[i]==y)
	return(ind)
}
