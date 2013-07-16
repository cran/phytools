# some utility functions
# written by Liam J. Revell 2011, 2012, 2013

# function rotates a node or multiple nodes
# written by Liam J. Revell 2013
rotateNodes<-function(tree,nodes,polytom=c(1,2),...){
	n<-length(tree$tip.label)
	if(nodes[1]=="all") nodes<-1:tree$Nnode+n
	for(i in 1:length(nodes)) tree<-rotate(tree,nodes[i],polytom)
	if(hasArg(reversible)) reversible<-list(...)$reversible
	else reversible<-TRUE
	if(reversible){ 
		ii<-which(tree$edge[,2]<=n)
		jj<-tree$edge[ii,2]
		tree$edge[ii,2]<-1:n
		tree$tip.label<-tree$tip.label[jj]
	}
	return(tree)
}

# function simulates random sampling from xbar, xvar, with sample sizes n
# written by Liam J. Revell 2012
sampleFrom<-function(xbar=0,xvar=1,n=1,randn=NULL,type="norm"){
	if(length(xvar)==1&&length(xbar)!=length(xvar)) xvar<-rep(xvar,length(xbar))
	if(!is.null(randn))
		for(i in 1:length(xbar)) n[i]<-floor(runif(n=1,min=randn[1],max=(randn[2]+1)))
	x<-vector()
	for(i in 1:length(xbar)){
		y<-rnorm(n=n[i],mean=xbar[i],sd=sqrt(xvar[i]))
   		names(y)<-rep(names(xbar)[i],length(y))
   		x<-c(x,y)
	}
	return(x)
}

# function adds a new tip to the tree
# written by Liam J. Revell 2012, 2013
bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL,position=0){
	if(is.null(where)) where<-length(tree$tip)+1
	if(where<=length(tree$tip.label)&&position==0){
		pp<-1e-12
		if(tree$edge.length[which(tree$edge[,2]==where)]<=1e-12){
			tree$edge.length[which(tree$edge[,2]==where)]<-2e-12
			ff<-TRUE
		} else ff<-FALSE
	} else pp<-position
	if(is.null(edge.length)&&is.ultrametric(tree)){
		H<-nodeHeights(tree)
		if(where==(length(tree$tip)+1)) edge.length<-max(H)
		else edge.length<-max(H)-H[tree$edge[,2]==where,2]+position
	}
	tip<-list(edge=matrix(c(2,1),1,2),
		tip.label=tip.label,
		edge.length=edge.length,
		Nnode=1)
		class(tip)<-"phylo"
	obj<-bind.tree(tree,tip,where=where,position=pp)
	if(where<=length(tree$tip.label)&&position==0){
		nn<-obj$edge[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label)),1]
		obj$edge.length[which(obj$edge[,2]==nn)]<-obj$edge.length[which(obj$edge[,2]==nn)]+1e-12
		obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label))]<-0
		obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tree$tip.label[where]))]<-0
	}
	return(obj)
}

# function collapses the subtree descended from node to a star tree
# written by Liam J. Revell 2013
collapse.to.star<-function(tree,node){
	tt<-splitTree(tree,split=list(node=node,bp=tree$edge.length[which(tree$edge[,2]==node)]))
	ss<-starTree(species=tt[[2]]$tip.label,branch.lengths=diag(vcv(tt[[2]])))
	ss$root.edge<-0
	tree<-paste.tree(tt[[1]],ss)
	return(tree)
}

# function returns the MRCA, or its height above the root, for a set of taxa (in tips)
# written by Liam Revell 2012, 2013
findMRCA<-function(tree,tips=NULL,type=c("node","height")){
	type<-type[1]
	if(is.null(tips)){ 
		X<-mrca(tree)
		if(type=="height"){
			H<-nodeHeights(tree)
			X<-apply(X,c(1,2),function(x,y,z) y[which(z==x)[1]],y=H,z=tree$edge)
		}
		return(X)
	} else {
		H<-nodeHeights(tree)
		X<-sapply(tips,function(x,y,z) sapply(y,fastMRCA,sp1=x,tree=z),y=tips,z=tree)
		Y<-apply(X,c(1,2),function(x,y,z) y[which(z==x)[1]],y=H,z=tree$edge)
		if(type=="height") return(Y[which.min(Y)]) else return(X[which.min(Y)])
	}
}

# fast pairwise MRCA function
# written by Liam Revell 2012
fastMRCA<-function(tree,sp1,sp2){
	x<-match(sp1,tree$tip.label)
	y<-match(sp2,tree$tip.label)
	a<-Ancestors(tree,x)
	b<-Ancestors(tree,y)
	z<-a%in%b
	return(a[min(which(z))])
}

# function reorders simmap tree
# written Liam Revell 2011, 2013
reorderSimmap<-function(tree,order="cladewise"){
	x<-reorder(tree,order)
	o<-whichorder(x$edge[,2],tree$edge[,2])
	x$mapped.edge<-tree$mapped.edge[o,]
	x$maps<-tree$maps[o]
	return(x)
}

# function whichorder
# written by Liam Revell 2011, 2013
whichorder<-function(x,y) sapply(x,function(x,y) which(x==y),y=y)

# function works like extract.clade in ape but will preserve a discrete character mapping
# written by Liam J. Revell 2013
extract.clade.simmap<-function(tree,node){
	x<-getDescendants(tree,node)
	x<-x[x<=length(tree$tip.label)]
	drop.tip.simmap(tree,tree$tip.label[-x])
}

# function gets all subtrees that cannot be further subdivided into two clades of >= clade.size
# written by Liam J. Revell 2013
getCladesofSize<-function(tree,clade.size=2){
	n<-length(tree$tip.label)
	nn<-1:(tree$Nnode+n)
	ndesc<-function(tree,node){
		x<-getDescendants(tree,node)
		sum(x<=length(tree$tip.label))
	}
	dd<-setNames(sapply(nn,ndesc,tree=tree),nn)
	aa<-n+1 # root
	nodes<-vector()
	while(length(aa)){
		bb<-lapply(aa,function(x,tree) tree$edge[which(tree$edge[,1]==x),2],tree=tree)
		cc<-lapply(bb,function(x) dd[as.character(x)])
		gg<-sapply(cc,function(x,cs) any(x<cs),cs=clade.size)
		nodes<-c(nodes,aa[gg])
		aa<-unlist(bb[!gg])
	}
	trees<-lapply(nodes,extract.clade,phy=tree)
	class(trees)<-"multiPhylo"
	return(trees)
}

# function to merge mapped states
# written by Liam J. Revell 2013
mergeMappedStates<-function(tree,old.states,new.state){
	if(class(tree)=="multiPhylo"){
		tree<-unclass(tree)
		tree<-lapply(tree,mergeMappedStates,old.states=old.states,new.state=new.state)
		class(tree)<-"multiPhylo"
	} else {
		maps<-tree$maps
		rr<-function(map,oo,nn){ 
			for(i in 1:length(map)) if(names(map)[i]%in%oo) names(map)[i]<-nn
			map
		}
		mm<-function(map){
			if(length(map)>1){
				new.map<-vector()
				j<-1
				new.map[j]<-map[1]
				names(new.map)[j]<-names(map)[1]
				for(i in 2:length(map)){
					if(names(map)[i]==names(map)[i-1]){ 
						new.map[j]<-map[i]+new.map[j]
						names(new.map)[j]<-names(map)[i]
					} else {
						j<-j+1
						new.map[j]<-map[i]
						names(new.map)[j]<-names(map)[i]
					}
				}
				map<-new.map
			}
			map
		}
		maps<-lapply(maps,rr,oo=old.states,nn=new.state)
		maps<-lapply(maps,mm)
		mapped.edge<-tree$mapped.edge
		mapped.edge<-cbind(rowSums(mapped.edge[,colnames(mapped.edge)%in%old.states]),
			mapped.edge[,setdiff(colnames(mapped.edge),old.states)])
		colnames(mapped.edge)<-c(new.state,setdiff(colnames(tree$mapped.edge),old.states))
		tree$maps<-maps
		tree$mapped.edge<-mapped.edge
	}
	return(tree)
}

# function to rescale simmap style trees
# written by Liam J. Revell 2012, 2013
rescaleSimmap<-function(tree,totalDepth=1.0){
	if(class(tree)=="multiPhylo"){
		trees<-unclass(tree)
		trees<-lapply(trees,rescaleSimmap,totalDepth)
		class(trees)<-"multiPhylo"
		return(trees)
	} else if(class(tree)=="phylo"){
		h<-max(nodeHeights(tree))
		s<-totalDepth/h
		tree$edge.length<-tree$edge.length*s
		maps<-lapply(tree$maps,"*",s)
		tree$maps<-maps
		tree$mapped.edge<-tree$mapped.edge*s
		return(tree)
	} else message("tree should be an object of class \"phylo\" or \"multiPhylo\"")
}

# function to summarize the results of stochastic mapping
# written by Liam J. Revell 2013
describe.simmap<-function(tree,...){
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(check.equal)) check.equal<-list(...)$check.equal
	else check.equal<-FALSE
	if(hasArg(message)) message<-list(...)$message
	else message<-TRUE
	if(class(tree)=="multiPhylo"){
		if(check.equal){
			TT<-sapply(tree,function(x,y) sapply(y,all.equal.phylo,x),y=tree)
			check<-all(TT)
			if(!check) warning("some trees not equal")
		}
		YY<-getStates(tree)
		states<-sort(unique(as.vector(YY)))
		ZZ<-t(apply(YY,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=length(tree)))
		XX<-countSimmap(tree,states,FALSE)
		XX<-XX[,-(which(as.vector(diag(-1,length(states)))==-1)+1)]
		AA<-t(sapply(unclass(tree),function(x) c(colSums(x$mapped.edge),sum(x$edge.length))))
		colnames(AA)[ncol(AA)]<-"total"
		BB<-getStates(tree,type="tips")
		CC<-t(apply(BB,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=length(tree)))
		if(message){
			cat(paste(length(tree),"trees with a mapped discrete character with states:\n",paste(states,collapse=", "),"\n\n"))
			cat(paste("trees have",colMeans(XX)["N"],"changes between states on average\n\n"))
			cat(paste("changes are of the following types:\n"))
			aa<-t(as.matrix(colMeans(XX)[2:ncol(XX)])); rownames(aa)<-"x->y"; print(aa)
			cat(paste("\nmean total time spent in each state is:\n"))
			print(matrix(c(colMeans(AA),colMeans(AA[,1:ncol(AA)]/AA[,ncol(AA)])),2,ncol(AA),byrow=TRUE,dimnames=list(c("raw","prop"),c(colnames(AA)))))
			cat("\n")
		}
		if(plot){
			plot(tree[[1]],edge.width=2,no.margin=TRUE,label.offset=0.015*max(nodeHeights(tree[[1]])),...)
			nodelabels(pie=ZZ,piecol=1:length(states),cex=0.6)
			tips<-CC
			tiplabels(pie=tips,piecol=1:length(states),cex=0.5)
			L<-list(count=XX,times=AA,tips=tips,ace=ZZ,legend=setNames(states,palette()[1:length(states)]))
			if(message) invisible(L) else return(L)
		} else {
			L<-list(count=XX,times=AA,ace=ZZ)
			if(message) invisible(L) else return(L)
		}
	} else if(class(tree)=="phylo"){
		XX<-countSimmap(tree,message=FALSE)
		YY<-getStates(tree)
		states<-sort(unique(YY))
		AA<-setNames(c(colSums(tree$mapped.edge),sum(tree$edge.length)),c(colnames(tree$mapped.edge),"total"))
		AA<-rbind(AA,AA/AA[length(AA)]); rownames(AA)<-c("raw","prop")
		if(message){
			cat(paste("1 tree with a mapped discrete character with states:\n",paste(states,collapse=", "),"\n\n"))
			cat(paste("tree has",XX$N,"changes between states\n\n"))
			cat(paste("changes are of the following types:\n"))
			print(XX$Tr)
			cat(paste("\nmean total time spent in each state is:\n"))
			print(AA)
			cat("\n")
		}
		if(plot){
			plotSimmap(tree,colors=setNames(palette()[1:length(states)],states),pts=FALSE)
			L<-list(N=XX$N,Tr=XX$Tr,times=AA,states=YY,legend=setNames(states,palette()[1:length(states)]))
			if(message) invisible(L) else return(L)
		} else { 
			L<-list(N=XX$N,Tr=XX$Tr,times=AA,states=YY)
			if(message) invisible(L) else return(L)
		}
	}
}

# function to get states at internal nodes from simmap style trees
# written by Liam J. Revell 2013
getStates<-function(tree,type=c("nodes","tips")){
	type<-type[1]
	if(class(tree)=="multiPhylo"){
		tree<-unclass(tree)
		y<-sapply(tree,getStates,type=type)
	} else if(class(tree)=="phylo"){ 
		if(type=="nodes"){
			y<-setNames(sapply(tree$maps,function(x) names(x)[1]),tree$edge[,1])
			y<-y[as.character(length(tree$tip)+1:tree$Nnode)]
		} else if(type=="tips"){
			y<-setNames(sapply(tree$maps,function(x) names(x)[length(x)]),tree$edge[,2])
			y<-setNames(y[as.character(1:length(tree$tip.label))],tree$tip.label)
		}
	} else stop("tree should be an object of class 'phylo' or 'multiPhylo'")
	return(y)
}

# function counts transitions from a mapped history
# written by Liam J. Revell 2013
countSimmap<-function(tree,states=NULL,message=TRUE){
	if(class(tree)=="multiPhylo"){
		ff<-function(zz){
 			XX<-countSimmap(zz,states,message)
			setNames(c(XX$N,as.vector(t(XX$Tr))),c("N",
			sapply(rownames(XX$Tr),paste,colnames(XX$Tr),sep=",")))
		}
		tree<-unclass(tree)
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

# function to re-root a phylogeny along an edge
# written by Liam Revell 2011-2013
reroot<-function(tree,node.number,position){
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	tt<-splitTree(tree,list(node=node.number,bp=position))
	p<-tt[[1]]; d<-tt[[2]]
	p<-root(p,outgroup="NA",resolve.root=T)
	bb<-which(p$tip.label=="NA")
	ee<-p$edge.length[which(p$edge[,2]==bb)]
	p$edge.length[which(p$edge[,2]==bb)]<-0
	cc<-p$edge[which(p$edge[,2]==bb),1]
	dd<-setdiff(p$edge[which(p$edge[,1]==cc),2],bb)
	p$edge.length[which(p$edge[,2]==dd)]<-p$edge.length[which(p$edge[,2]==dd)]+ee
	tt<-paste.tree(p,d)
	return(tt)
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

