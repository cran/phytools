# some utility functions
# written by Liam J. Revell 2011, 2012, 2013

# function counts transitions from a mapped history
# written by Liam J. Revell 2013
countSimmap<-function(tree,states=NULL,message=TRUE){
	if(class(tree)=="multiPhylo"){
		ff<-function(zz){
 			XX<-countSimmap(zz,states,message)
			setNames(c(XX$N,as.vector(t(XX$Tr))),c("N",
			sapply(rownames(XX$Tr),paste,colnames(XX$Tr),sep=",")))
		}
		XX<-t(sapply(tree,ff))	
		if(!message) return(XX)
		else return(list(Tr=XX,message=
			c("Column N is the total number of character changes on the tree",
			"Other columns give transitions x,y from x->y")))
	} else if(class(tree)=="phylo") {
		n<-sum(sapply(tree$maps,length))-nrow(tree$edge)
		if(is.null(states)) states<-colnames(tree$mapped.edge)
		m<-length(states)	
		TT<-matrix(NA,m,m,dimnames=list(states,states))
		gg<-function(map,a,b){
			if(length(map)==1) zz<-0
			else {
				zz<-0; i<-2
				while(i<=length(map)){
					if(names(map)[i]==b&&names(map)[i-1]==a) zz<-zz+1
				i<-i+1
				}
			} 
			return(zz)
		}
		for(i in 1:m) for(j in 1:m)
			if(i==j) TT[i,j]<-0
			else TT[i,j]<-sum(sapply(tree$maps,gg,a=states[i],b=states[j]))
		if(!message) return(list(N=n,Tr=TT))
		else return(list(N=n,Tr=TT,message=c(
			"N is the total number of character changes on the tree",
			"Tr gives the number of transitions from row state->column state")))
	}
}

# function returns random state with probability given by y
# written by Liam J. Revell 2013 (replaces earlier version)
rstate<-function(y){
	if(length(y)==1) return(names(y)[1])
	else return(names(which(rmultinom(1,1,y/sum(y))[,1]==1)))
}

# function to match nodes between trees
# written by Liam J. Revell 2012, 2013
matchNodes<-function(tr1,tr2,method=c("descendants","distances"),...){
	method<-method[1]
	method<-matchType(method,c("descendants","distances"))
	if(method=="descendants"){
		desc.tr1<-lapply(1:tr1$Nnode+length(tr1$tip),function(x) extract.clade(tr1,x)$tip.label)
		names(desc.tr1)<-1:tr1$Nnode+length(tr1$tip)
		desc.tr2<-lapply(1:tr2$Nnode+length(tr2$tip),function(x) extract.clade(tr2,x)$tip.label)
		names(desc.tr2)<-1:tr2$Nnode+length(tr2$tip)
		Nodes<-matrix(NA,length(desc.tr1),2,dimnames=list(NULL,c("tr1","tr2")))
		for(i in 1:length(desc.tr1)){
			Nodes[i,1]<-as.numeric(names(desc.tr1)[i])
			for(j in 1:length(desc.tr2))
				if(all(desc.tr1[[i]]%in%desc.tr2[[j]])&&all(desc.tr2[[j]]%in%desc.tr1[[i]]))
					Nodes[i,2]<-as.numeric(names(desc.tr2)[j])
		}
	} else if(method=="distances"){
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-1e-6
		if(hasArg(corr)) corr<-list(...)$corr
		else corr<-FALSE
		if(corr) tr1$edge.length<-tr1$edge.length/max(nodeHeights(tr1))
		if(corr) tr2$edge.length<-tr2$edge.length/max(nodeHeights(tr2))
		D1<-dist.nodes(tr1)[1:length(tr1$tip),1:tr1$Nnode+length(tr1$tip)]
		D2<-dist.nodes(tr2)[1:length(tr2$tip),1:tr2$Nnode+length(tr2$tip)]
		rownames(D1)<-tr1$tip.label
		rownames(D2)<-tr2$tip.label
		common.tips<-intersect(tr1$tip.label,tr2$tip.label)
		D1<-D1[common.tips,]
		D2<-D2[common.tips,]
		Nodes<-matrix(NA,tr1$Nnode,2,dimnames=list(NULL,c("tr1","tr2")))
		for(i in 1:tr1$Nnode){
			if(corr) z<-apply(D2,2,function(X,y) cor(X,y),y=D1[,i])
			else z<-apply(D2,2,function(X,y) 1-sum(abs(X-y)),y=D1[,i])
			Nodes[i,1]<-as.numeric(colnames(D1)[i])
			if(any(z>=(1-tol))){
				a<-as.numeric(names(which(z>=(1-tol))))
				if(length(a)==1) Nodes[i,2]<-a
				else {
					Nodes[i,2]<-a[1]
					warning("polytomy detected; some node matches may be arbitrary")
				}
			}
		}
	}
	return(Nodes)
}

# function rounds the branch lengths of the tree & applies rounding to simmap tree
# written by Liam J. Revell 2012
roundBranches<-function(tree,digits=0){
	if(class(tree)=="multiPhylo"){
		trees<-lapply(tree,roundBranches)
		class(trees)<-"multiPhylo"
		return(trees)
	} else {
		tree$edge.length<-round(tree$edge.length,digits)
		if(!is.null(tree$maps)){
			for(i in 1:nrow(tree$edge)){
				temp<-tree$maps[[i]]/sum(tree$maps[[i]])
				tree$maps[[i]]<-temp*tree$edge.length[i]
			}
		}
		if(!is.null(tree$mapped.edge)){
			a<-vector()
			for(i in 1:nrow(tree$edge)) a<-c(a,names(tree$maps[[i]]))
			a<-unique(a)
			tree$mapped.edge<-matrix(data=0,length(tree$edge.length),length(a),dimnames=list(apply(tree$edge,1,function(x) paste(x,collapse=",")),state=a))
			for(i in 1:length(tree$maps)) for(j in 1:length(tree$maps[[i]])) tree$mapped.edge[i,names(tree$maps[[i]])[j]]<-tree$mapped.edge[i,names(tree$maps[[i]])[j]]+tree$maps[[i]][j]
		}
		return(tree)
	}
}

# function applies the branch lengths of a reference tree to a second tree, including mappings
# written by Liam J. Revell 2012
applyBranchLengths<-function(tree,edge.length){
	if(class(tree)=="multiPhylo"){
		trees<-lapply(tree,applyBranchLengths,edge.length=edge.length)
		class(trees)<-"multiPhylo"
		return(trees)
	} else {
		tree$edge.length<-edge.length
		if(!is.null(tree$maps)){
			for(i in 1:nrow(tree$edge)){
				temp<-tree$maps[[i]]/sum(tree$maps[[i]])
				tree$maps[[i]]<-temp*tree$edge.length[i]
			}
		}
		if(!is.null(tree$mapped.edge)){
			a<-vector()
			for(i in 1:nrow(tree$edge)) a<-c(a,names(tree$maps[[i]]))
			a<-unique(a)
			tree$mapped.edge<-matrix(data=0,length(tree$edge.length),length(a),dimnames=list(apply(tree$edge,1,function(x) paste(x,collapse=",")),state=a))
			for(i in 1:length(tree$maps)) for(j in 1:length(tree$maps[[i]])) tree$mapped.edge[i,names(tree$maps[[i]])[j]]<-tree$mapped.edge[i,names(tree$maps[[i]])[j]]+tree$maps[[i]][j]
		}
		return(tree)
	}
}

# function to compute phylogenetic VCV using joint Pagel's lambda
# written by Liam Revell 2011
phyl.vcv<-function(X,C,lambda){
	C<-lambda.transform(lambda,C)
	invC<-solve(C)
	a<-matrix(colSums(invC%*%X)/sum(invC),ncol(X),1)
	A<-matrix(rep(a,nrow(X)),nrow(X),ncol(X),byrow=T)
	V<-t(X-A)%*%invC%*%(X-A)/(nrow(C)-1)
	return(list(C=C,R=V,alpha=a))
}

# lambda transformation of C
# written by Liam Revell 2011
lambda.transform<-function(lambda,C){
	if(lambda==1) return(C)
	else {
		V<-diag(diag(C))
		C<-C-V
		C.lambda<-(V+lambda*C)
		return(C.lambda)
	}
}

# likelihood function for joint estimation of lambda for multiple traits
# written by Liam Revell 2011/2012
likMlambda<-function(lambda,X,C){
	# compute R, conditioned on lambda
	temp<-phyl.vcv(X,C,lambda);
	C<-temp$C; R<-temp$R; a<-temp$alpha
	# prep
	n<-nrow(X); m<-ncol(X); D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	y<-as.matrix(as.vector(X))
	# compute the log-likelihood
	kronRC<-kronecker(R,C)
	logL<--t(y-D%*%a)%*%solve(kronRC,y-D%*%a)/2-n*m*log(2*pi)/2-determinant(kronRC,logarithm=TRUE)$modulus/2
	return(logL)
}

# function matches data to tree
# written by Liam J. Revell 2011
matchDatatoTree<-function(tree,x,name){
	if(is.matrix(x)) x<-x[,1]
	if(is.null(names(x))){
		if(length(x)==length(tree$tip)){
			print(paste(name,"has no names; assuming x is in the same order as tree$tip.label"))
			names(x)<-tree$tip.label
		} else
			stop(paste(name,"has no names and is a different length than tree$tip.label"))
	}
	if(any(is.na(match(names(x),tree$tip.label)))){
		print(paste("some species in",name,"are missing from tree, dropping missing taxa from",name))
		x<-x[intersect(tree$tip.label,names(x))]
	}
	return(x)
}

# function matches tree to data
# written by Liam J. Revell 2011
matchTreetoData<-function(tree,x,name){
	if(any(is.na(match(tree$tip.label,names(x))))){
		print(paste("some species in tree are missing from",name,", dropping missing taxa from the tree"))
		tree<-drop.tip(tree,setdiff(tree$tip.label,names(x)))
	}
	if(any(is.na(x))){
		print(paste("some data in",name,"given as 'NA', dropping corresponding species from tree"))
		tree<-drop.tip(tree,names(which(is.na(x))))
	}
	return(tree)
}

# function finds the maximum value of Pagel's lambda
# written by Liam J. Revell 2011
maxLambda<-function(tree){
	if(is.ultrametric(tree)){
		H<-nodeHeights(tree)
		return(max(H[,2])/max(H[,1]))
	} else return(1)
}

# function reorders the columns of mapped.edge from a set of simmap trees
# written by Liam J. Revell 2013
orderMappedEdge<-function(trees,ordering=NULL){
	f1<-function(tree,ordering){
		mapped.edge<-matrix(0,nrow(tree$mapped.edge),length(ordering),
			dimnames=list(rownames(tree$mapped.edge),ordering))
		mapped.edge[,colnames(tree$mapped.edge)]<-tree$mapped.edge
		tree$mapped.edge<-mapped.edge
		return(tree)
	}
	f2<-function(tree) colnames(tree$mapped.edge)
	if(class(trees)=="phylo") states<-colnames(trees$mapped.edge)
	else if(class(trees)=="multiPhylo") states<-unique(as.vector(sapply(trees,f2)))
	else stop("trees should be an object of class 'phylo' or 'multiPhylo'")
	if(length(ordering)>1) if(length(intersect(states,ordering))<length(states)){
		warning("not all states represented in input ordering; setting to default")
		ordering<-NULL
	}
	if(is.null(ordering)) ordering<-"alphabetical"
	if(length(ordering)==1){
		ordering<-matchType(ordering,c("alphabetical","numerical"))
		if(ordering=="alphabetical") ordering<-sort(states)
		else if(ordering=="numerical") ordering<-as.character(sort(as.numeric(states)))
	}
	if(class(trees)=="phylo") trees<-f1(trees,ordering)
	else { 
		trees<-lapply(trees,f1,ordering=ordering)
		class(trees)<-"multiPhylo"
	}
	return(trees)
}

# function gets sister node numbers or names
# written by Liam J. Revell 2013
getSisters<-function(tree,node,mode=c("number","label")){
	mode<-mode[1]
	if(is.character(node)) node<-match(node,c(tree$tip.label,tree$node.label))
	sisters<-tree$edge[which(tree$edge[,1]==tree$edge[which(tree$edge[,2]==node),1]),2]
	sisters<-setdiff(sisters,node)
	if(mode=="number") return(sisters)
	else if(mode=="label"){
		result<-list()
		n<-length(tree$tip.label)
		if(is.null(tree$node.label)&&any(sisters>n)) result$nodes<-sisters[which(sisters>n)] 
		else if(any(sisters>n)) result$nodes<-tree$node.label[sisters[which(sisters>n)]-n]
		if(any(sisters<=n)) result$tips<-tree$tip.label[sisters[which(sisters<=n)]]
		return(result)
	}
}

# gets descendant node numbers
# written by Liam Revell 2012, 2013
getDescendants<-function(tree,node,curr=NULL){
	if(is.null(curr)) curr<-vector()
	daughters<-tree$edge[which(tree$edge[,1]==node),2]
	curr<-c(curr,daughters)
	if(length(curr)==0&&node<=length(tree$tip.label)) curr<-node
	w<-which(daughters>=length(tree$tip))
	if(length(w)>0) for(i in 1:length(w)) 
		curr<-getDescendants(tree,daughters[w[i]],curr)
	return(curr)
}

# function computes vcv for each state, and stores in array
# written by Liam J. Revell 2011/2012
multiC<-function(tree){
	n<-length(tree$tip); m<-ncol(tree$mapped.edge)
	# compute separate C for each state
	mC<-list()
	for(i in 1:m){
		mtree<-list(edge=tree$edge,Nnode=tree$Nnode,tip.label=tree$tip.label,edge.length=tree$mapped.edge[,i])
		class(mtree)<-"phylo"
		mC[[i]]<-vcv.phylo(mtree)
	}
	return(mC)
}

# function pastes subtree onto tip
# written by Liam Revell 2011
paste.tree<-function(tr1,tr2){
	if(length(tr2$tip)>1){ 
		temp<-tr2$root.edge; tr2$root.edge<-NULL
		tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]<-tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]+temp
	}
	tr.bound<-bind.tree(tr1,tr2,where=which(tr1$tip.label=="NA"))	
	return(tr.bound)
}

# function drops entire clade
# written by Liam Revell 2011
drop.clade<-function(tree,tip){
	tree<-drop.tip(tree,tip,trim.internal=FALSE)
	while(sum(tree$tip.label=="NA")>1){
		tree<-drop.tip(tree,"NA",trim.internal=FALSE)
	}
	return(tree)
}

# match type
# written by Liam J. Revell 2012
matchType<-function(type,types){
	for(i in 1:length(types))
		if(all(strsplit(type,split="")[[1]]==strsplit(types[i],split="")[[1]][1:length(strsplit(type,split="")[[1]])]))
			type=types[i]
	return(type)
}

# wraps around MatrixExp
# written by Liam Revell 2011
expm<-function(Y){
	Z<-MatrixExp(Y); dimnames(Z)<-dimnames(Y)
	return(Z)
}
	
# function reorders simmap tree
# written Liam Revell 2011
reorderSimmap<-function(tree,order="cladewise"){
	ntree<-reorder(tree,order)
	o<-whichorder(ntree$edge[,2],tree$edge[,2])
	ntree$mapped.edge<-tree$mapped.edge[o,]
	ntree$maps<-tree$maps[o]
	return(ntree)
}

# function 'untangles' (or attempts to untangle) a tree with crossing branches
# written by Liam J. Revell 2013
untangle<-function(tree,method=c("reorder","read.tree")){
	method<-method[1]
	if(!is.null(tree$maps)) simmap<-TRUE
	else simmap<-FALSE
	if(method=="reorder"){
		if(simmap) tree<-reorderSimmap(reorderSimmap(tree,"pruningwise"))
		else tree<-reorder(reorder(tree,"pruningwise"))
	} else if(method=="read.tree"){
		if(simmap){
			stop("Option 'read.tree' does not presently work for SIMMAP style trees")
			# tree<-read.simmap(text=write.simmap(tree))
		} else tree<-read.tree(text=write.tree(tree))
	}
	return(tree)
}

