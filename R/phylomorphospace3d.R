# phylomorphospace3d: projection of a tree into three dimensional morphospace
# written by Liam J. Revell 2012

phylomorphospace3d<-function(tree,X,A=NULL,label=TRUE,control=list()){
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	# control list
	con=list(spin=TRUE,axes=TRUE,box=TRUE,simple.axes=FALSE,lwd=1,ftype="reg")
	con[(namc<-names(control))]<-control
	if(con$simple.axes) con$box<-con$axes<-FALSE
	con$ftype<-which(c("off","reg","b","i","bi")==con$ftype)-1
	if(is.null(A)) A<-apply(X,2,function(x,tree) anc.ML(tree,x)$ace,tree=tree)
	else A<-A[as.character(1:tree$Nnode+length(tree$tip)),]
	x<-y<-z<-matrix(NA,nrow(tree$edge),2)
	X<-X[tree$tip.label,]
	for(i in 1:length(tree$tip)){
		x[tree$edge[,2]==i,2]<-X[i,1]
		y[tree$edge[,2]==i,2]<-X[i,2]
		z[tree$edge[,2]==i,2]<-X[i,3]
	}
	for(i in length(tree$tip)+1:tree$Nnode){
		x[tree$edge[,1]==i,1]<-x[tree$edge[,2]==i,2]<-A[as.character(i),1]
		y[tree$edge[,1]==i,1]<-y[tree$edge[,2]==i,2]<-A[as.character(i),2]
		z[tree$edge[,1]==i,1]<-z[tree$edge[,2]==i,2]<-A[as.character(i),3]
	}
	if(is.null(colnames(X))) colnames(X)<-c("x","y","z")
	plot3d(rbind(X,A),xlab=colnames(X)[1],ylab=colnames(X)[2],zlab=colnames(X)[3],axes=con$axes,box=con$box)
	spheres3d(X,radius=0.02*mean(apply(X,2,range)))
	for(i in 1:nrow(tree$edge)) segments3d(x[i,],y[i,],z[i,],lwd=con$lwd)
	ms<-colMeans(X)
	rs<-apply(rbind(X,A),2,range)
	if(con$simple.axes){
		lines3d(x=rs[,1],y=c(rs[1,2],rs[1,2]),z=c(rs[1,3],rs[1,3]))
		lines3d(x=c(rs[1,1],rs[1,1]),y=rs[,2],z=c(rs[1,3],rs[1,3]))
		lines3d(x=c(rs[1,1],rs[1,1]),y=c(rs[1,2],rs[1,2]),z=rs[,3])
	}
	rs<-rs[2,]-rs[1,]
	for(i in 1:length(tree$tip)){
		adj<-0.03*rs*(2*(X[i,]>ms)-1)
		if(con$ftype) text3d(X[i,]+adj,texts=tree$tip.label[i],font=con$ftype)
	}
	if(con$spin) play3d(spin3d(axis=c(0,0,1),rpm=10),duration=5)
}	
