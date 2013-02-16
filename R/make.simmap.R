# function creates a stochastic character mapped tree as a modified "phylo" object
# written by Liam Revell 2011, 2013

make.simmap<-function(tree,x,model="SYM",nsim=1){
	# reorder to cladewise
	tree<-reorder.phylo(tree,"cladewise")
	# fit discrete character model
	XX<-ace(x,tree,type="discrete",model=model)
	# get Q matrix
	Q<-matrix(XX$rates[XX$index.matrix],nrow(XX$index.matrix),ncol(XX$index.matrix),dimnames=list(colnames(XX$lik.anc),colnames(XX$lik.anc)))
	L<-XX$lik.anc; rownames(L)<-length(tree$tip)+1:nrow(L)
	diag(Q)<--colSums(Q,na.rm=TRUE)
	# create "multiPhylo" object
	mtrees<-list(); class(mtrees)<-"multiPhylo"
	for(i in 1:nsim){
		# create the map tree object
		mtree<-tree; mtree$maps<-list()
		# now we want to simulate the node states on the tree
		node.states<-matrix(NA,nrow(tree$edge),ncol(tree$edge))
		node.states[which(tree$edge[,1]==(length(tree$tip)+1)),1]<-rstate(L[1,])
		for(j in 1:nrow(tree$edge)){
			if(tree$edge[j,2]>length(tree$tip)){
				p<-expm(Q*tree$edge.length[j])[node.states[j,1],]*L[as.character(tree$edge[j,2]),]
				node.states[j,2]<-rstate(p/max(p))
				node.states[which(tree$edge[,1]==tree$edge[j,2]),1]<-node.states[j,2]
			} else {
				node.states[j,2]<-x[tree$tip.label[tree$edge[j,2]]]
			}
			# now simulate on the branches
			if(tree$edge.length[j]==0){ 
				map<-vector(); map[1]<-tree$edge.length[j]; names(map)[1]<-node.states[j,2]
			} else {
				accept=FALSE
				while(!accept){
					time=0; state<-node.states[j,1]; new.state<-state; dt<-vector(); map<-vector(); k<-1
					while(time<tree$edge.length[j]){
						dt[1]<-time
						dt[2]<-dt[1]+rexp(n=1,rate=-Q[state,state])
						if(dt[2]<tree$edge.length[j]) new.state<-rstate(Q[,state][-match(state,rownames(Q))]/sum(Q[,state][-match(state,rownames(Q))]))
						dt[2]<-min(dt[2],tree$edge.length[j])
						map[k]<-dt[2]-dt[1]; names(map)[k]<-state; k<-k+1
						state<-new.state; time<-dt[2]
					}
					if(names(map)[length(map)]==node.states[j,2]) accept=TRUE
				}
			}	
			mtree$maps[[j]]<-map
		}
		# now construct the matrix "mapped.edge" (for backward compatibility
		allstates<-vector()
		for(j in 1:nrow(mtree$edge)) allstates<-c(allstates,names(mtree$maps[[j]]))
		allstates<-unique(allstates)
		mtree$mapped.edge<-matrix(data=0,length(mtree$edge.length),length(allstates),dimnames=list(apply(mtree$edge,1,function(x) paste(x,collapse=",")),state=allstates))
		for(j in 1:length(mtree$maps)) for(k in 1:length(mtree$maps[[j]])) mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
		mtrees[[i]]<-mtree	
	}
	if(nsim==1) mtrees<-mtrees[[1]]
	return(mtrees)
}

