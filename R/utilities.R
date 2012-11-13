# some utility functions used by phyl.cca and phyl.pca
# written by Liam J. Revell 2011/2012

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

# function
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

# function
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

# function
# written by Liam J. Revell 2011

maxLambda<-function(tree){
	if(is.ultrametric(tree)){
		H<-nodeHeights(tree)
		return(max(H[,2])/max(H[,1]))
	} else return(1)
}

# function to match nodes between trees
# written by Liam J. Revell 2012

matchNodes<-function(tr1,tr2){
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
	return(Nodes)
}

