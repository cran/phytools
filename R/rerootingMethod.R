# function to compute the marginal posterior probabilities for nodes using the rerooting method
# written by Liam J. Revell 2013

rerootingMethod<-function(tree,x,model=c("ER","SYM")){
	model<-model[1]
	n<-length(tree$tip.label)
	XX<-matrix(NA,tree$Nnode,length(unique(x)))
	rownames(XX)<-1:tree$Nnode+n
	colnames(XX)<-sort(unique(x))
	YY<-apeAce(tree,x,model=model)
	XX[1,]<-YY$lik.anc[1,]
	Q<-matrix(YY$rates[YY$index.matrix],ncol(XX),ncol(XX),dimnames=list(colnames(XX),colnames(XX)))
	diag(Q)<--colSums(Q,na.rm=TRUE)
	for(i in 2:tree$Nnode){
		tt<-reroot(tree,node.number=i+n,position=tree$edge.length[which(tree$edge[,2]==(i+n))])
		XX[i,]<-apeAce(tt,x,model=model)$lik.anc[1,]
	}
	return(list(loglik=YY$loglik,Q=Q,marginal.anc=XX))
}
	
