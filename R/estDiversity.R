# function computes an estimate of the standing diversity in each category given by x at each node
# written by Liam J. Revell 2011, 2013

estDiversity<-function(tree,x){
	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	# first get the node heights
	tree<-reorder(tree,"cladewise")
	node.heights<-nodeHeights(tree)
	# now reconstruct x on tree
	anc<-ace(x,tree,type="discrete")$lik.anc
	rownames(anc)<-1:tree$Nnode+length(tree$tip)
	# now loop through every node above the root
	D<-matrix(0,nrow(anc),ncol(anc),dimnames=dimnames(anc))
	D[1,]<-rep(0,ncol(anc))
	message("Please wait. . . . Warning - this may take a while!")
	for(i in 2:nrow(anc)){
		t<-node.heights[match(rownames(anc)[i],tree$edge[,2]),2]
		ind<-node.heights[,1]<t&node.heights[,2]>t
		edges<-matrix(tree$edge[ind,],length(tree$edge[ind,])/2,2)
		heights<-matrix(node.heights[ind,],nrow(edges),ncol(edges))
		for(j in 1:nrow(edges)){
			tr<-reroot(tree,edges[j,2],t-heights[j,1])
			D[i,]<-D[i,]+ace(x,tr,type="discrete")$lik.anc[1,]
		}
		D[i,]<-D[i,]*anc[i,]
		if(i%%10==0) message(paste("Completed",i,"nodes"))
	}
	d<-rowSums(D)
}

