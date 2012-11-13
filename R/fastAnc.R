# function does fast estimation of ML ancestral states using ace
# written by Liam J. Revell 2012

fastAnc<-function(tree,x){
	if(!is.binary.tree(tree)) btree<-multi2di(tree)
	else btree<-tree
	M<-btree$Nnode
	N<-length(btree$tip)
	anc<-vector()
	for(i in 1:M+N){
   		anc[i-N]<-ace(x,multi2di(root(btree,node=i)),method="pic")$ace[1]
   		names(anc)[i-N]<-i
 	}
	if(!is.binary.tree(tree)){
		ancNames<-matchNodes(tree,btree)
		anc<-anc[as.character(ancNames[,2])]
		names(anc)<-ancNames[,1]
	}
	return(anc)
}
