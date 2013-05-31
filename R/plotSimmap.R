# function plots a stochastic character mapped tree
# written by Liam Revell 2011, 2013

plotSimmap<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,pts=TRUE,node.numbers=FALSE,mar=NULL,add=FALSE,offset=NULL,direction="rightwards",type="phylogram"){
	if(class(tree)=="multiPhylo"){
		par(ask=TRUE)
		for(i in 1:length(tree)) plotSimmap(tree[[i]],colors=colors,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers,mar,add,offset,direction,type)
	} else {
		# check font
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!ftype) fsize=0 
		# check colors
		if(is.null(colors)){
			st<-sort(unique(unlist(sapply(tree$maps,names))))
			colors<-palette()[1:length(st)]
			names(colors)<-st
			if(length(st)>1){
				cat("no colors provided. using the following legend:\n")
				print(colors)
			}
		}
		# check tree
		if(class(tree)!="phylo") stop("tree should be object of class \"phylo\"")
		if(is.null(tree$maps)) stop("tree should contain mapped states on edges.")
		# swap out "_" character for spaces (assumes _ is a place holder)
		tree$tip.label<-gsub("_"," ",tree$tip.label)
		# get margin
		if(is.null(mar)) mar=rep(0.1,4)
		if(type=="phylogram"){
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
			par(mar=mar)
			if(!add) plot.new()
			if(fsize*max(strwidth(cw$tip.label))<1.0){
				c<-(1-fsize*max(strwidth(cw$tip.label)))/max(node.height)
				cw$edge.length<-c*cw$edge.length
				cw$maps<-lapply(cw$maps,function(x) x<-c*x)
				node.height<-c*node.height
			} else message("Font size too large to properly rescale tree to window.")
			if(!add){
				if(direction=="leftwards") plot.window(xlim=c(max(node.height)+fsize*max(strwidth(cw$tip.label)),0),ylim=c(1,max(Y)))
				else plot.window(xlim=c(0,max(node.height)+fsize*max(strwidth(cw$tip.label))),ylim=c(1,max(Y)))
			}
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
				symbols(0,mean(Y[cw$edge[cw$edge[,1]==(length(cw$tip)+1),2]]),rectangles=matrix(c(1.2*fsize*strwidth(as.character(length(cw$tip)+1)),1.4*fsize*strheight(as.character(length(cw$tip)+1))),1,2),inches=FALSE,bg="white",add=TRUE)
				text(0,mean(Y[cw$edge[cw$edge[,1]==(length(cw$tip)+1),2]]),length(cw$tip)+1,cex=fsize)
				for(i in 1:nrow(cw$edge)){
					x<-node.height[i,2]
					if(cw$edge[i,2]>length(tree$tip)){
						symbols(x,Y[cw$edge[i,2]],rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,bg="white",add=TRUE)
						text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
					}
				}
			}
			if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
			pos<-if(direction=="leftwards") 2 else 4
			for(i in 1:n) if(ftype) text(node.height[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,offset=offset,cex=fsize,font=ftype)
			# reset margin
			par(mar=c(5,4,4,2)+0.1)
		} else if(type=="fan"){
			plotFan(tree,colors,fsize,ftype,lwd,mar,add)
		}
	}
}	

# function to plot simmap tree in type "fan"
# written by Liam J. Revell 2013
plotFan<-function(tree,colors,fsize,ftype,lwd,mar,add){
	cat("\nNote: type=\"fan\" is in development.\nMany options of type=\"phylogram\" are not yet available.\n\n")
	# reorder
	cw<-reorder(tree)
	pw<-reorder(tree,"pruningwise")
	# count nodes and tips
	n<-length(cw$tip)
	m<-cw$Nnode 
	# get Y coordinates on uncurved space
	Y<-vector(length=m+n)
	Y[cw$edge[cw$edge[,2]<=length(cw$tip),2]]<-1:n
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
		Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
	}
	Y<-setNames(Y/max(Y)*2*pi,1:(n+m))
	Y<-cbind(Y[as.character(tree$edge[,2])],Y[as.character(tree$edge[,2])])
	R<-nodeHeights(cw)
	# now put into a circular coordinate system
	x<-R*cos(Y)
	y<-R*sin(Y)
	# optimize x & y limits
	par(mar=mar)
	offsetFudge<-1.37 # empirically determined
	offset<-0
	pp<-par("pin")[1]
 	sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+offsetFudge*offset*fsize*strwidth("W",units="inches") 
	alp<-optimize(function(a,H,sw,pp) (2*a*1.04*max(H)+2*sw-pp)^2,H=R,sw=sw,pp=pp,interval=c(0,1e6))$minimum 
	xylim<-c(-max(R)-sw/alp,max(R)+sw/alp)
	# plot tree
	if(!add) plot.new()
	plot.window(xlim=xylim,ylim=xylim,asp=1,)
	# plot radial lines (edges)
	for(i in 1:nrow(cw$edge)){
		maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
		xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
		yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
		for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],lwd=lwd,lend=2)
	}
	# plot circular lines
	for(i in 1:m+n){
		r<-R[match(i,cw$edge)]
		a1<-min(Y[which(cw$edge==i)])
		a2<-max(Y[which(cw$edge==i)])
		draw.arc(0,0,r,a1,a2,lwd=lwd,col=colors[names(cw$maps[[match(i,cw$edge[,1])]])[1]])
	}
	# plot labels
	for(i in 1:n){
		ii<-which(cw$edge[,2]==i)
		aa<-Y[ii,2]/(2*pi)*360
		adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
		tt<-if(aa>90&&aa<270) paste(cw$tip.label[i]," ",sep="") else paste(" ",cw$tip.label[i],sep="")
		aa<-if(aa>90&&aa<270) 180+aa else aa
		if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
	}
}
