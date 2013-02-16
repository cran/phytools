# function performs ancestral character estimation under the threshold model
# written by Liam J. Revell 2012,2013

ancThresh<-function(tree,x,ngen=1000,sequence=NULL,method="mcmc",control=list(),...){
	
	# check method
	if(method!="mcmc") stop(paste(c("do not recognize method =",method,",quitting")))

	# check x
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)){
		X<-x[tree$tip.label,]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, column order in x")
 			seq<-colnames(X)
		} else seq<-sequence
	} else if(is.vector(x)){
		x<-x[tree$tip.label]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, using alphabetical or numerical order")
 			seq<-sort(levels(as.factor(x)))
		} else seq<-sequence
		X<-to.matrix(x,seq)
	}
	# row scale X
	X<-X/apply(X,1,sum)
	X<-X[,seq] # order columns by seq

	# ok, now set starting thresholds
	th<-c(1:length(seq))-1; names(th)<-seq
	x<-to.vector(X)
	l<-sapply(x,function(x) runif(n=1,min=th[x]-1,max=th[x])) # set plausible starting liability

	# for MCMC
	n<-length(tree$tip)
	m<-length(th)
	npar<-tree$Nnode+n+m-2

	# populate control list
	PrA<-matrix(1/m,tree$Nnode,m,dimnames=list(1:tree$Nnode+n,seq))
	if(!is.null(control$pr.anc)){
		if(!is.matrix(control$pr.anc)){
			message("**** NOTE: prior on ancestral states must be in matrix form; using default prior")
			control$pr.anc<-NULL
		} else {
			control$pr.anc<-control$pr.anc[,seq,drop=FALSE]
			PrA[rownames(control$pr.anc),]<-control$pr.anc
			control$pr.anc<-PrA	
		}
	}
	con=list(sample=1000,
		propliab=0.5*max(nodeHeights(tree)),
		propthresh=0.05*max(nodeHeights(tree)),
		pr.anc=PrA,
		pr.th=0.01,
		burnin=round(0.2*ngen),
		plot=TRUE,
		print=TRUE,
		piecol=NULL,
		tipcol="input")
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]

	# now set ancestral liabilities, by first picking ancestral states from their prior
	temp<-apply(con$pr.anc,1,rstate)
	# assign random liabilities consistent with the starting thresholds
	a<-sapply(temp,function(x) runif(n=1,min=th[x]-1,max=th[x]))

	# now change the upper limit of th to Inf
	th[length(th)]<-Inf

	# compute some matrices & values
	V<-vcvPhylo(tree)
	# check to make sure that V will be non-singular
	if(any(tree$edge.length<=(10*.Machine$double.eps)))
		stop("some branch lengths are 0 or nearly zero")
	invV<-solve(V)
	detV<-determinant(V,logarithm=TRUE)$modulus[1]
	lik1<-likLiab(l,a,V,invV,detV)+log(probMatch(X,l,th,seq))
	
	# store
	A<-matrix(NA,ngen/con$sample+1,tree$Nnode,dimnames=list(NULL,n+1:tree$Nnode))
	B<-matrix(NA,ngen/con$sample+1,m+2,dimnames=list(NULL,c("gen",names(th),"logLik")))
	C<-matrix(NA,ngen/con$sample+1,tree$Nnode+n,dimnames=list(NULL,c(tree$tip.label,1:tree$Nnode+n)))
	A[1,]<-sapply(a,threshState,thresholds=th)
	B[1,]<-c(0,th,lik1)
	C[1,]<-c(l[tree$tip.label],a[as.character(1:tree$Nnode+n)])

	# run MCMC
	message("MCMC starting....")
	for(i in 1:ngen){
		lik1<-likLiab(l,a,V,invV,detV)+log(probMatch(X,l,th,seq))
		d<-i%%npar
		if(ngen>=1000) if(i%%1000==0) if(con$print) message(paste("gen",i))
		ap<-a; lp<-l; thp<-th
		if(d<=tree$Nnode&&d!=0){
			# update node liabilities
			ind<-d%%tree$Nnode; if(ind==0) ind<-tree$Nnode
			ap[ind]<-a[ind]+rnorm(n=1,sd=sqrt(con$propliab))
		} else {
			if((d>tree$Nnode&&d<=(tree$Nnode+n))||(npar==(tree$Nnode+n)&&d==0)){
				# update tip liabilities
				if(d==0) ind<-n
				else { ind<-(d-tree$Nnode)%%n; if(ind==0) ind<-n }
				lp[ind]<-l[ind]+rnorm(n=1,sd=sqrt(con$propliab))
			} else {
				# update thresholds
				if(d) ind<-(d-tree$Nnode-n)%%m+1
				else ind<-m-1
				thp[ind]<-bounce(th[ind],rnorm(n=1,sd=sqrt(con$propthresh)),c(th[ind-1],th[ind+1]))
			}
		}
		lik2<-likLiab(lp,ap,V,invV,detV)+log(probMatch(X,lp,thp,seq))
		p.odds<-min(c(1,exp(lik2+logPrior(sapply(ap,threshState,thresholds=thp),thp,con)-lik1-logPrior(sapply(a,threshState,thresholds=th),th,con))))

		if(p.odds>runif(n=1)){
			a<-ap; l<-lp; th<-thp
			logL<-lik2
		} else logL<-lik1
		if(i%%con$sample==0){ 
			A[i/con$sample+1,]<-sapply(a,threshState,thresholds=th)
			B[i/con$sample+1,]<-c(i,th[colnames(B)[1+1:m]],logL)
			C[i/con$sample+1,]<-c(l[tree$tip.label],a[as.character(1:tree$Nnode+n)])
		}
	}
	mcmc<-as.data.frame(A)
	param<-as.data.frame(B)
	liab<-as.data.frame(C)
	ace<-matrix(0,tree$Nnode,m,dimnames=list(colnames(A),seq))
	burnin<-which(param[,"gen"]==con$burnin)
	for(i in 1:tree$Nnode){
		temp<-summary(mcmc[burnin:nrow(mcmc),i])/(nrow(mcmc)-burnin+1)
		ace[i,names(temp)]<-temp
	}
	if(con$plot) plotThresh(tree,X,list(ace=ace,mcmc=mcmc,par=param,liab=liab),burnin=con$burnin,piecol=con$piecol,tipcol=con$tipcol,...)
	return(list(ace=ace,mcmc=mcmc,par=param,liab=liab))
}

# plots ancestral states from the threshold model
# written by Liam J. Revell 2012

plotThresh<-function(tree,x,mcmc,burnin=NULL,piecol,tipcol="input",...){

	# plot tree
	plot(tree,no.margin=TRUE,edge.width=2,...)
	
	# pull matrices from mcmc
	ace<-mcmc$ace
	liab<-mcmc$liab
	param<-mcmc$par

	# get burnin
	if(is.null(burnin)) burnin<-round(0.2*max(param[,"gen"]))
	burnin<-which(param[,"gen"]==burnin)

	# check x
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)) X<-x[tree$tip.label,]
	else if(is.vector(x)){
		x<-x[tree$tip.label]
		X<-to.matrix(x,names(piecol))
	}
	# row scale X
	X/apply(X,1,sum)->X

	# plot node labels
	nodelabels(pie=ace,piecol=piecol[colnames(ace)],cex=0.6)

	# plot tip labels
	if(tipcol=="input") tiplabels(pie=X,piecol=piecol[colnames(X)],cex=0.6)
	else if(tipcol=="estimated") {
		XX<-matrix(NA,nrow(liab),length(tree$tip),dimnames=list(rownames(liab),colnames(liab)[1:length(tree$tip)]))
		for(i in 1:nrow(liab)) XX[i,]<-sapply(liab[i,1:length(tree$tip)],threshState,thresholds=param[i,1:ncol(X)+1])
		X<-t(apply(XX,2,function(x) summary(factor(x,levels=colnames(X)))))
		tiplabels(pie=X/rowSums(X),piecol=piecol[colnames(X)],cex=0.6)
	}
}

# computes DIC for threshold model
# written by Liam J. Revell 2012

threshDIC<-function(tree,x,mcmc,burnin=NULL,sequence=NULL,method="pD"){

	# check x
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)){
		X<-x[tree$tip.label,]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, column order in x")
 			seq<-colnames(X)
		} else seq<-sequence
	} else if(is.vector(x)){
		x<-x[tree$tip.label]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, using alphabetical or numerical order")
 			seq<-sort(levels(as.factor(x)))
		} else seq<-sequence
		X<-to.matrix(x,seq)
	}
	# row scale X
	X<-X/apply(X,1,sum)
	X<-X[,seq] # order columns by seq

	# convert burnin to starting row
	if(is.null(burnin)) burnin<-0.2*max(mcmc$par[,"gen"])
	start<-which(mcmc$par[,"gen"]==burnin)+1

	# compute
	thBar<-colMeans(mcmc$par[start:nrow(mcmc$par),2:(ncol(mcmc$par)-1)])
	liabBar<-colMeans(mcmc$liab[start:nrow(mcmc$liab),])
	Dtheta<--2*(likLiab(liabBar[tree$tip.label],liabBar[as.character(1:tree$Nnode+length(tree$tip))],vcvPhylo(tree),solve(vcvPhylo(tree)),determinant(vcvPhylo(tree),logarithm=TRUE)$modulus[1])+log(probMatch(X[tree$tip.label,],liabBar[tree$tip.label],thBar,seq)))
	D<--2*mcmc$par[start:nrow(mcmc$par),"logLik"]
	Dbar<-mean(D)
	if(method=="pD"){		
		pD<-Dbar-Dtheta
		DIC<-pD+Dbar
		result<-setNames(c(Dbar,Dtheta,pD,DIC),c("Dbar","Dhat","pD","DIC"))
	} else if(method=="pV"){
		pV<-var(D)/2
		DIC<-pV+Dbar
		result<-setNames(c(Dbar,Dtheta,pV,DIC),c("Dbar","Dhat","pV","DIC"))
	}
	return(result)
}

# internal functions for ancThresh, plotThresh, and threshDIC

# returns a state based on position relative to thresholds
threshState<-function(x,thresholds){
	t<-c(-Inf,thresholds,Inf)
	names(t)[length(t)]<-names(t)[length(t)-1] 
	i<-1; while(x>t[i]) i<-i+1
	return(names(t)[i])
}

# likelihood function for the liabilities
likLiab<-function(l,a,V,invV,detV){
	y<-c(l,a[2:length(a)])-a[1]
	logL<--y%*%invV%*%y/2-nrow(V)*log(2*pi)/2-detV/2
	return(logL)
}

# function for the log-prior
logPrior<-function(a,t,control){
	pp<-sum(log(diag(control$pr.anc[names(a),a])))+
		if(length(t)>2) sum(dexp(t[2:(length(t)-1)],rate=control$pr.th,log=TRUE)) else 0				
	return(pp)		
}

# check if the liability predictions match observed data
allMatch<-function(x,l,thresholds){
	result<-all(sapply(l,threshState,thresholds=thresholds)==x)
	if(!is.na(result)) return(result)
	else return(FALSE)
}

# check if the liability predictions match observed data & return a probability
# (this allows states to be uncertain)
probMatch<-function(X,l,thresholds,sequence){
	Y<-to.matrix(sapply(l,threshState,thresholds=thresholds),sequence)
	return(prod(rowSums(X*Y)))
}

# bounds parameter by bouncing
bounce<-function(start,step,bounds){
	x<-start+step
	while(x>bounds[2]||x<bounds[1]){
		if(x>bounds[2]) x<-2*bounds[2]-x
		if(x<bounds[1]) x<-2*bounds[1]-x
	}
	return(x)
}

# convert vector of x to binary matrix
to.matrix<-function(x,seq){
	X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
	for(i in 1:length(seq)) X[x==seq[i],i]<-1
	return(X)
}

# convert binary matrix to vector
to.vector<-function(X) apply(X,1,rstate)
