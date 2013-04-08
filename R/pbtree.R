# function to simulate a pure-birth phylogenetic tree or trees
# written by Liam J. Revell 2011/2012

pbtree<-function(b=1,n=NULL,t=NULL,scale=NULL,nsim=1,ape=TRUE){
	if(nsim>1){
		trees<-replicate(nsim,pbtree(b,n,t,scale),simplify=FALSE)
		class(trees)<-"multiPhylo"
		return(trees)
	} else {
		if(!is.null(t))
			stop("time stop not yet implemented")
		else {
			node<-n+1
			edge<-matrix(c(node,NA,node,NA),2,2,byrow=T)
			edge.length<-c(0,0)
			node<-node+1; tip<-0
			while(nrow(edge)<(2*n-2)){
				# new edge
				o<-is.na(edge[,2])
				p<-which(o)
				q<-sample(p)[1]
				edge[q,2]<-node
				edge<-rbind(edge,matrix(c(node,NA,node,NA),2,2,byrow=T))
				node<-node+1
				l<-rexp(n=1,sum(o)*b)
				edge.length[p]<-edge.length[p]+l
				edge.length<-c(edge.length,rep(0,2))
			}
			o<-is.na(edge[,2])
			p<-which(o)
			l<-rexp(n=1,sum(o)*b)			
			edge.length[p]<-edge.length[p]+l		
			edge[is.na(edge[,2]),2]<-tip+1:sum(is.na(edge[,2]))
			tip.label<-paste("t",1:n,sep="")
			tree<-list(edge=edge,edge.length=edge.length,tip.label=tip.label,Nnode=n-1)
			class(tree)<-"phylo"
			if(!is.null(scale)){
				h<-max(nodeHeights(tree))
				tree$edge.length<-scale*tree$edge.length/h
			}
			if(ape) tree<-read.tree(text=write.tree(tree))
			return(tree)
		}
	}
}
