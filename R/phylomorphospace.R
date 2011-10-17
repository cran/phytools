# This function creates a phylomorsphace plot of Sidlauskas (2006)
# Written by Liam J. Revell 2010 (updated 3-8-2011)
phylomorphospace<-function(tree,X,A=NULL,label=TRUE,control=list()){

	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	if(nrow(X)!=length(tree$tip)) stop("X must contain the same number of rows as species in tree.")
	if(ncol(X)!=2) warning("X has more than 2 columns.  Using only the first 2 columns.")

	# check to see if A should be estimated
	if(is.null(A)){
		A<-matrix(NA,tree$Nnode,2,dimnames=list(as.character(length(tree$tip.label)+1:tree$Nnode),colnames(X)))
		for(i in 1:2) A[,i]<-ace(X[,i],tree)$ace
	}

	# control list
	con=list(col.edge=matrix(data="black",nrow(tree$edge),1,dimnames=list(as.character(tree$edge[,2]),"color")),col.node=matrix(data="black",nrow(tree$edge)+1,1,dimnames=list(as.character(c(tree$edge[1,1],tree$edge[,2])),"color")))
	con[(namc<-names(control))]<-control
	con$col.edge<-as.matrix(con$col.edge); con$col.node<-as.matrix(con$col.node)

	# plot root state
	plot(x=A[1,1],y=A[1,2],xlim=range(c(X[,1],A[,1])),ylim=range(c(X[,2],A[,2])),xlab=colnames(X)[1],ylab=colnames(X)[2])
	
	# put X in tree order
	X<-X[tree$tip.label,]
	rownames(X)<-1:length(tree$tip)
	
	# check for node labels; create if necessary
	attach(tree); if(!exists("node.label")) tree$node.label<-as.character(length(tree$tip.label)+1:tree$Nnode); detach(tree)
	tree$node.label<-as.character(tree$node.label)
	A<-A[tree$node.label,]
	rownames(A)<-length(tree$tip)+1:tree$Nnode

	# bind X & A
	X<-rbind(X,A)

	# loop across the branches of the tree
	for(i in 1:nrow(tree$edge)){
		lines(x=c(X[as.character(tree$edge[i,1]),1],X[as.character(tree$edge[i,2]),1]),y=c(X[as.character(tree$edge[i,1]),2],X[as.character(tree$edge[i,2]),2]),col=con$col.edge[as.character(tree$edge[i,2]),1])
		points(x=X[as.character(tree$edge[i,1]),1],y=X[as.character(tree$edge[i,1]),2],col=con$col.node[as.character(tree$edge[i,1]),1],pch=16,cex=1.0)
		if(tree$edge[i,2]<=length(tree$tip.label)&&label)
			textxy(X=X[as.character(tree$edge[i,2]),1],Y=X[as.character(tree$edge[i,2]),2],labs=tree$tip.label[tree$edge[i,2]],cx=0.75)	
	}

	# plot larger points for the tips
	points(X[1:length(tree$tip),1],X[1:length(tree$tip),2],col=con$col.node[rownames(X)[1:length(tree$tip)],1],pch=16,cex=1.3)

}
