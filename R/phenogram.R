# function
# written by Liam J. Revell 2011

phenogram<-function(tree,x,fsize=1.0,ftype="reg"){
	# check tree
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo'")
	# check font
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	if(!ftype) fsize=0 
	H<-nodeHeights(tree)
	if(length(x)<(length(tree$tip)+tree$Nnode))
		x<-c(x,anc.ML(tree,x)$ace)
	X<-matrix(x[tree$edge],nrow(tree$edge),ncol(tree$edge))
	plot.new()
	plot.window(ylim=c(min(x),max(x)),xlim=c(min(H),max(H)+fsize*max(strwidth(tree$tip.label))))
	for(i in 1:nrow(H)){ 
	lines(H[i,],X[i,])
	if(tree$edge[i,2]<=length(tree$tip))
		if(fsize) text(tree$tip.label[tree$edge[i,2]],x=H[i,2]+0.02*max(H),y=X[i,2],cex=fsize,font=ftype)
	}
	axis(1); axis(2)
}

