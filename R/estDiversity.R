# function computes an estimate of the standing diversity in each category given by x at each node
# written by Liam J. Revell 2011

estDiversity<-function(tree,x){
	# check for & load "ape"	
	if(!require(ape)) stop("must first install 'ape' package.") # require ape	
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

# reroot tree along an edge
# written by Liam J. Revell 2011

reroot<-function(tree,node.number,position){
	if(node.number>length(tree$tip)){
		# first, re-root the tree at node.number
		tr<-root(tree,node=node.number,resolve.root=T)
		# now re-allocate branch length to the two edges descending from the new root node
		b<-sum(tr$edge.length[tr$edge==(length(tree$tip)+1)])
		tr$edge.length[tr$edge==(length(tree$tip)+1)]<-c(position,b-position)
	} else {
		# first, root the tree at the parent of node.number
		tr1<-root(tree,node=tree$edge[match(node.number,tree$edge[,2]),1])
		# now drop tip
		tr1<-drop.tip(tr1,tree$tip.label[node.number])
		# create phylo object
		tr2<-list(edge=matrix(c(3L,1L,3L,2L),2,2,byrow=T),tip.label=c(tree$tip.label[node.number],"NA"),edge.length=c(tree$edge.length[match(node.number,tree$edge[,2])]-position,position),Nnode=1)
		class(tr2)<-"phylo"
		tr<-bind.tree(tr2,tr1,where=which(tr2$tip.label=="NA"))
	}
	return(tr)
}

