## code to place a missing extant taxon into a tree using ML on continuous data
## written by Liam J. Revell 2014

locate.yeti<-function(tree,X,...){
	if(hasArg(method)) method<-list(...)$method
	else method<-"heuristic"
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(constraint)){
		if(method=="exhaustive") constraint<-list(...)$constraint
		else{
			cat("constraint only works with method==\"exhaustive\"\n")
			constraint<-tree$edge[,2]
		}
	} else constraint<-tree$edge[,2]
	if(!is.matrix(X)) X<-as.matrix(X)
	tip<-setdiff(rownames(X),tree$tip.label)
	if(ncol(X)>1){
		pca<-phyl.pca(tree,X[tree$tip.label,])
		obj<-phyl.vcv(X[tree$tip.label,],vcv(tree),1)
		X<-(X-matrix(rep(obj$a[,1],nrow(X)),nrow(X),ncol(X),byrow=TRUE))%*%pca$Evec
	}
	if(method=="heuristic"){
		trees<-list()
		for(i in 1:nrow(tree$edge)) trees[[i]]<-bind.tip(tree,tip,where=tree$edge[i,2],position=0.5*tree$edge.length[i])
		class(trees)<-"multiPhylo"
		lik.edge<-function(tree,XX){
			obj<-phyl.vcv(as.matrix(XX[tree$tip.label,]),vcv(tree),1)
			ll<-vector()
			for(i in 1:ncol(XX)) ll[i]<-sum(dmnorm(XX[tree$tip.label,i],mean=rep(obj$a[i,1],nrow(XX)),obj$C*obj$R[i,i],log=TRUE))
			sum(ll)
		}
		logL<-sapply(trees,lik.edge,XX=X)
		if(plot){
			ll<-logL
			ll[ll<=sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]]<-sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]
			layout(matrix(c(1,2),2,1),heights=c(0.95,0.05))
			plotBranchbyTrait(tree,ll,mode="edges",title="log(L)",show.tip.label=FALSE)
			plot.new()
			text(paste("Note: logL <=",round(min(ll),2),"set to",round(min(ll),2),"for visualization only"),x=0.5,y=0.5)
		}
		edge<-tree$edge[which(logL==max(logL)),2]
	}
	lik.tree<-function(position,tip,tree,edge,XX){
		tree<-bind.tip(tree,tip,where=edge,position=position)
		obj<-phyl.vcv(as.matrix(XX[tree$tip.label,]),vcv(tree),1)
		ll<-vector()
		for(i in 1:ncol(XX)) ll[i]<-sum(dmnorm(XX[tree$tip.label,i],mean=rep(obj$a[i,1],nrow(XX)),obj$C*obj$R[i,i],log=TRUE))
		sum(ll)
	}
	if(method=="heuristic"){
		ee<-edge
		if(tree$edge[which(tree$edge[,2]==edge),1]!=(length(tree$tip.label)+1)) ee<-c(ee,getAncestors(tree,node=edge,type="parent"))
		if(edge>length(tree$tip.label)) ee<-c(ee,tree$edge[which(tree$edge[,1]==edge),2])
	} else if(method=="exhaustive") ee<-c(1:length(tree$tip.label),2:tree$Nnode+length(tree$tip.label))
	ee<-intersect(ee,constraint)
	fit<-vector(mode="list",length=length(ee))
	for(i in 1:length(ee)) fit[[i]]<-optimize(lik.tree,interval=c(0,tree$edge.length[which(tree$edge[,2]==ee[i])]),tip=tip,tree=tree,edge=ee[i],XX=X,maximum=TRUE)
	logL<-sapply(fit,function(x) x$objective)
	if(method=="exhaustive"&&plot){
		ll<-sapply(fit,function(x) x$objective)
		ll[ll<=sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]]<-sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]
		layout(matrix(c(1,2),2,1),heights=c(0.95,0.05))
		plotBranchbyTrait(tree,ll,mode="edges",title="log(L)",show.tip.label=FALSE)
		plot.new()
		text(paste("Note: logL <=",round(min(ll),2),"set to",round(min(ll),2),"for visualization only"),x=0.5,y=0.5)
	}
	fit<-fit[[which(logL==max(logL))]]
	edge<-ee[which(logL==max(logL))]
	mltree<-bind.tip(tree,tip,where=edge,position=fit$maximum)
	mltree$logL<-fit$objective
	mltree
}

