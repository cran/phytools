# function to compute the marginal posterior probabilities for nodes using the rerooting method
# written by Liam J. Revell 2013

rerootingMethod<-function(tree,x,model=c("ER","SYM"),...){
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-NULL
	if(!is.matrix(model)) model<-model[1]
	n<-length(tree$tip.label)
	# if vector convert to binary matrix
	if(!is.matrix(x)){ 
		yy<-to.matrix(x,sort(unique(x)))
		if(is.null(tips)) tips<-FALSE
	} else { 
		if(is.null(tips)) tips<-TRUE
		yy<-x
	}
	yy<-yy[tree$tip.label,]
	yy<-yy/rowSums(yy)
	XX<-matrix(NA,tree$Nnode+n,ncol(yy))
	rownames(XX)<-1:(tree$Nnode+n)
	colnames(XX)<-colnames(yy)
	YY<-apeAce(tree,yy,model=model)
	XX[n+1,]<-YY$lik.anc[1,]
	Q<-matrix(c(0,YY$rates)[YY$index.matrix+1],ncol(XX),ncol(XX),dimnames=list(colnames(XX),colnames(XX)))
	diag(Q)<--colSums(Q,na.rm=TRUE)
	for(i in 1:(tree$Nnode+n)){
		if(i!=(n+1)){
			if(i>n||tips){
				tt<-reroot(tree,node.number=i,position=tree$edge.length[which(tree$edge[,2]==i)])
				XX[i,]<-apeAce(tt,yy,model=model,fixedQ=Q)$lik.anc[1,]
			}
		}
	}
	rownames(XX)[1:n]<-tree$tip.label
	XX<-if(tips) XX else XX[1:tree$Nnode+n,]
	return(list(loglik=YY$loglik,Q=Q,marginal.anc=XX))
}
