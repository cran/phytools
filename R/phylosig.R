# function for computing phylogenetic signal by the lambda (Pagel 1999) of K (Blomberg et al. 2003) methods
# written by Liam J. Revell 2011

phylosig<-function(tree,x,method="K",test=FALSE,nsim=1000){
	# some minor error checking
	if(class(tree)!="phylo") stop("tree object must be of class 'phylo.'")
	if(is.matrix(x)) x<-x[,1]
	if(is.null(names(x))){
		if(length(x)==length(tree$tip)){
			print("x has no names; assuming x is in the same order as tree$tip.label")
			names(x)<-tree$tip.label
		} else
			stop("x has no names and is a different length than tree$tip.label")
	}
	if(any(is.na(match(tree$tip.label,names(x))))){
		print("some species in tree are missing from data, dropping missing taxa from the tree")
		tree<-drop.tip(tree,tree$tip.label[-match(names(x),tree$tip.label)])
	}
	if(any(is.na(match(names(x),tree$tip.label)))){
		print("some species in data are missing from tree, dropping missing taxa from the data")
		x<-x[tree$tip.label]
	}
	if(any(is.na(x))){
		print("some data given as 'NA', dropping corresponding species from tree")
		tree<-drop.tip(tree,names(which(is.na(x))))
	}
	# done error handling
	if(method=="K"){
		C<-vcv.phylo(tree)
		x<-x[rownames(C)]
		invC<-solve(C)
		n<-nrow(C)
		a<-sum(invC%*%x)/sum(invC)
		K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invC%*%(x-a)))/((sum(diag(C))-n/sum(invC))/(n-1)) # calculate K
		if(!test)
			return(as.numeric(K)) # return K
		else {
			P=0.0
			simX<-x
			for(i in 1:nsim){
				a<-sum(invC%*%simX)/sum(invC)
				simK<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invC%*%(simX-a)))/((sum(diag(C))-n/sum(invC))/(n-1))
				if(simK>=K) P<-P+1/nsim # calculate P-value for randomization test
				simX<-sample(simX) # randomize x
			}
			return(list(K=as.numeric(K),P=P)) # return K & P
		}
	} else if(method=="lambda"){
		# function to compute C with lambda
		lambda.transform<-function(C,lambda){
			dC<-diag(diag(C))
			C<-lambda*(C-dC)+dC
			return(C)
		}
		# likelihood function
		likelihood<-function(theta,C,y){
			Cl<-lambda.transform(C,theta)
			invCl<-solve(Cl)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			one<-matrix(1,n,1)
			a<-as.numeric(sum(invCl%*%y)/sum(invCl))
			sig2<-as.numeric(t(y-a)%*%invCl%*%(y-a)/n)
			logL<--t(y-a)%*%(1/sig2*invCl)%*%(y-a)/2-n*log(2*pi)/2-determinant(sig2*Cl,logarithm=TRUE)$modulus/2
			return(logL)
		}
		C<-vcv.phylo(tree)
		maxLambda<-max(C)/max(C[upper.tri(C)])
		res<-optimize(f=likelihood,interval=c(0,maxLambda),y=x,C=C,maximum=TRUE) # optimize lambda
		if(!test)
			return(list(lambda=res$maximum,logL=res$objective[1,1])) # return lambda and log-likelihood
		else {
			logL0<-likelihood(theta=0,C=C,y=x) # compute likelihood of lambda=0
			P<-1-as.numeric(pchisq(2*(res$objective[1,1]-logL0),df=1)) # P-value
			return(list(lambda=res$maximum,logL=res$objective[1,1],P=P)) # return lambda, logL, and P-value
		}
	} else
		stop(paste("do not recognize method = \"",method,"\"; methods are \"K\" and \"lambda\"",sep=""))
}
