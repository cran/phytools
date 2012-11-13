# function performs ancestral character estimation under the threshold model
# written by Liam J. Revell 2012

ancThresh<-function(tree,x,ngen=1000,sequence=NULL,method="mcmc",control=list()){

	# check method
	if(method!="mcmc") stop(paste(c("do not recognize method =",method,",quitting")))

	# likelihood function for the liabilities
	lik<-function(l,a,V,invV,detV){
		y<-c(l,a[2:length(a)])-a[1]
		logL<--y%*%invV%*%y/2-nrow(V)*log(2*pi)/2-detV/2
		return(logL)
	}

	# function for the log-prior (presently returns 0, i.e. flat prior)
	logPrior<-function(a,t){
		pp<-sum(log(diag(con$pr.anc[names(a),a])))+
			if(length(t)>2) sum(dexp(t[2:(length(t)-1)],rate=con$pr.th,log=TRUE)) else 0				
		return(pp)		
	}

	# returns a state based on position relative to thresholds
	threshState<-function(x,thresholds){
		t<-c(-Inf,thresholds,Inf)
		names(t)[length(t)]<-names(t)[length(t)-1] 
		i<-1; while(x>t[i]) i<-i+1
		return(names(t)[i])
	}

	# check if the liabilities match the predicted
	allMatch<-function(x,l,thresholds){
		result<-all(sapply(l,threshState,thresholds=thresholds)==x)
		if(!is.na(result)) return(result)
		else return(FALSE)
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

	# get the sequence, or use alphanumeric order
	if(is.null(sequence)){
		message("**** NOTE: no sequence provided, using alphabetical or numerical order")
 		seq<-sort(levels(as.factor(x)))
	} else seq<-sequence

	# ok, now pull out the unique states in x & set starting thresholds
	x<-x[tree$tip.label]
	th<-c(1:length(seq))-1; names(th)<-seq
	if(is.factor(x)){ temp<-names(x); as.character(x)->x; names(x)<-temp } 	
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
		piecol=NULL)
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
	invV<-solve(V)
	detV<-determinant(V,logarithm=TRUE)$modulus[1]
	lik1<-lik(l,a,V,invV,detV)+log(allMatch(x,l,th))
	
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
		lik1<-lik(l,a,V,invV,detV)+log(allMatch(x,l,th))
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
		lik2<-lik(lp,ap,V,invV,detV)+log(allMatch(x,lp,thp))
		p.odds<-min(c(1,exp(lik2+logPrior(sapply(ap,threshState,thresholds=thp),thp)-lik1-logPrior(sapply(a,threshState,thresholds=th),th))))

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
	par<-as.data.frame(B)
	liab<-as.data.frame(C)
	ace<-matrix(0,tree$Nnode,m,dimnames=list(colnames(A),seq))
	for(i in 1:tree$Nnode){
		burnin<-which(par[,"gen"]==con$burnin)
		temp<-summary(mcmc[burnin:nrow(mcmc),i])/(nrow(mcmc)-burnin+1)
		ace[i,names(temp)]<-temp
	}
	if(con$plot){
		plot(tree,no.margin=TRUE,show.tip.label=FALSE)
		if(is.null(con$piecol)){ 
			con$piecol<-1:length(th)
			con$piecol<-con$piecol[seq]
		}
		nodelabels(pie=ace,piecol=con$piecol,cex=0.6)
		P<-matrix(0,n,m,dimnames=list(tree$tip.label,names(th)))
		for(i in 1:n) P[i,which(colnames(P)==x[i])]<-1
		tiplabels(pie=P,piecol=con$piecol[seq],cex=0.6)
	}
	return(list(ace=ace,mcmc=mcmc,par=par,liab=liab))
}

# computes DIC for threshold model
# written by Liam J. Revell 2012

threshDIC<-function(tree,x,mcmc,burnin){

	# likelihood function for the liabilities
	lik<-function(l,a,V,invV,detV){
		y<-c(l,a[2:length(a)])-a[1]
		logL<--y%*%invV%*%y/2-nrow(V)*log(2*pi)/2-detV/2
		return(logL)
	}

	# returns a state based on position relative to thresholds
	threshState<-function(x,thresholds){
		t<-c(-Inf,thresholds,Inf)
		names(t)[length(t)]<-names(t)[length(t)-1] 
		i<-1; while(x>t[i]) i<-i+1
		return(names(t)[i])
	}

	# check if the liabilities match the predicted
	allMatch<-function(x,l,thresholds){
		result<-all(sapply(l,threshState,thresholds=thresholds)==x)
		if(!is.na(result)) return(result)
		else return(FALSE)
	}

	# convert burnin to starting row
	start<-which(mcmc$par[,"gen"]==burnin)+1

	# compute
	Dbar<-mean(mcmc$par[start:nrow(mcmc$par),"logLik"])
	thBar<-colMeans(mcmc$par[start:nrow(mcmc$par),2:(ncol(mcmc$par)-1)])
	liabBar<-colMeans(mcmc$liab[start:nrow(mcmc$liab),])
	Dtheta<-lik(liabBar[tree$tip.label],liabBar[as.character(1:tree$Nnode+length(tree$tip))],vcvPhylo(tree),solve(vcvPhylo(tree)),determinant(vcvPhylo(tree),logarithm=TRUE)$modulus[1])+
			log(allMatch(x[tree$tip.label],liabBar[tree$tip.label],thBar))
	pD<-Dbar-Dtheta
	DIC<-pD+Dbar
	return(as.numeric(DIC))

}
