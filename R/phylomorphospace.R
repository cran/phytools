# this funciton creates a phylomorphospace plot (Sidlauskas 2006)
# written by Liam J. Revell 2010-13

phylomorphospace<-function(tree,X,A=NULL,label=TRUE,control=list(),...){

	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	if(nrow(X)!=length(tree$tip)) stop("X must contain the same number of rows as species in tree.")
	if(is.null(rownames(X))){
		warning("X is missing row names; assuming order of tip labels.")
		rownames(X)<-tree$tip.label
	}
	if(ncol(X)!=2){ 
		warning("X has more than 2 columns.  Using only the first 2 columns.")
		X<-X[,1:2]
	}

	# get ancestral states
	if(is.null(A)) A<-apply(X,2,fastAnc,tree=tree)

	# control list
	con=list(col.edge=setNames(rep("black",nrow(tree$edge)),as.character(tree$edge[,2])),
		col.node=setNames(rep("black",max(tree$edge)),as.character(1:max(tree$edge))))
	con[(namc<-names(control))]<-control

	# get optional argument
	if(hasArg(node.by.map)) node.by.map<-list(...)$node.by.map
	else node.by.map<-FALSE

	# set xlim & ylim
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-range(c(X[,1],A[,1]))
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-range(c(X[,2],A[,2]))

	# set xlab & ylab
	if(hasArg(xlab)) xlab<-list(...)$xlab
	else xlab<-colnames(X)[1]
	if(hasArg(ylab)) ylab<-list(...)$ylab
	else ylab<-colnames(X)[2]

	# set font size for tip labels
	if(hasArg(fsize)) fsize<-0.75*list(...)$fsize
	else fsize<-0.75
	
	# check if colors for stochastic mapping
	if(hasArg(colors)) colors<-list(...)$colors
	else if(!is.null(tree$maps)) colors<-setNames(palette()[1:ncol(tree$mapped.edge)],sort(colnames(tree$mapped.edge)))

	# set lwd
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-if(is.null(tree$maps)) 1 else 2

	# do some bookkeeping
	aa<-setNames(c(X[tree$tip.label,1],A[,1]),c(1:length(tree$tip.label),rownames(A)))
	bb<-setNames(c(X[tree$tip.label,2],A[,2]),c(1:length(tree$tip.label),rownames(A)))
	XX<-matrix(aa[as.character(tree$edge)],nrow(tree$edge),2)
	YY<-matrix(bb[as.character(tree$edge)],nrow(tree$edge),2)

	# plot projection
	plot(x=A[1,1],y=A[1,2],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,pch=16,cex=0.1,col="white")
	if(is.null(tree$maps)){
		for(i in 1:nrow(XX)) lines(XX[i,],YY[i,],col=con$col.edge[as.character(tree$edge[i,2])],lwd=lwd)
	} else {
		for(i in 1:nrow(XX)){
			xx<-tree$maps[[i]]/sum(tree$maps[[i]])*(XX[i,2]-XX[i,1])
			yy<-tree$maps[[i]]/sum(tree$maps[[i]])*(YY[i,2]-YY[i,1])
			cc<-names(tree$maps[[i]])
			x<-XX[i,1]; y<-YY[i,1]
			for(j in 1:length(xx)){
				lines(c(x,x+xx[j]),c(y,y+yy[j]),col=colors[cc[j]],lwd=lwd)
				x<-x+xx[j]; y<-y+yy[j]
			}
		}
		if(node.by.map){
			zz<-c(getStates(tree,type="tips"),getStates(tree))
			names(zz)[1:length(tree$tip.label)]<-sapply(names(zz)[1:length(tree$tip.label)],function(x,y) which(y==x),y=tree$tip.label)
			con$col.node<-setNames(colors[zz],names(zz))
		}
	}
	zz<-c(tree$edge[1,1],tree$edge[tree$edge[,2]>length(tree$tip.label),2])
	points(c(XX[1,1],XX[tree$edge[,2]>length(tree$tip.label),2]),c(YY[1,1],YY[tree$edge[,2]>length(tree$tip.label),2]),pch=16,cex=1.0,col=con$col.node[as.character(zz)])
	zz<-tree$edge[tree$edge[,2]<=length(tree$tip.label),2]
	points(XX[tree$edge[,2]<=length(tree$tip.label),2],YY[tree$edge[,2]<=length(tree$tip.label),2],pch=16,cex=1.3,col=con$col.node[as.character(zz)])
	zz<-sapply(1:length(tree$tip.label),function(x,y) which(x==y),y=tree$edge[,2])
	if(label) textxy(XX[zz,2],YY[zz,2],labs=tree$tip.label,cx=fsize)
}
	
		
	

	
