# function creates a phenogram
# written by Liam J. Revell 2011

phenogram<-function(tree,x,fsize=1.0,ftype="reg",colors=NULL){
	# check tree
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	# check font
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	if(!ftype) fsize=0 
	H<-nodeHeights(tree)
	if(length(x)<(length(tree$tip)+tree$Nnode))
		x<-c(x,anc.ML(tree,x)$ace)
	else
		x<-c(x[tree$tip.label],x[as.character(length(tree$tip)+1:tree$Nnode)])
	x[1:length(tree$tip)]<-x[tree$tip.label]
	names(x)[1:length(tree$tip)]<-1:length(tree$tip)
	X<-matrix(x[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
	plot.new()
	plot.window(ylim=c(min(x),max(x)),xlim=c(min(H),max(H)+fsize*max(strwidth(tree$tip.label))))
	if(is.null(tree$maps)){
		for(i in 1:nrow(H)){ 
			lines(H[i,],X[i,])
			if(tree$edge[i,2]<=length(tree$tip))
				if(fsize) text(tree$tip.label[tree$edge[i,2]],x=H[i,2]+0.02*max(H),y=X[i,2],cex=fsize,font=ftype)
		}
	} else {
		if(is.null(colors)){ 
			colors<-palette()
			names(colors)<-as.character(1:8)
		}
		for(i in 1:nrow(H)){
			y<-H[i,1]
			m<-diff(X[i,])/diff(H[i,])
			for(j in 1:length(tree$maps[[i]])){
				a<-c(y,y+tree$maps[[i]][j])
				b<-m*(a-H[i,1])+X[i,1]
				lines(a,b,col=colors[names(tree$maps[[i]])[j]],lwd=2)
				y<-a[2]
			}
			if(tree$edge[i,2]<=length(tree$tip))
				if(fsize) text(tree$tip.label[tree$edge[i,2]],x=H[i,2]+0.02*max(H),y=X[i,2],cex=fsize,font=ftype)
		}
	}
	axis(1); axis(2); title(xlab="time",ylab="phenotype")
}

