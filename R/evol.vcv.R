# This function fits two or more different evolutionary variance-covariance matrices (rate matrices) 
# to different parts of a phylogenetic tree
# Takes a tree with a binary or multistate character "painted" on the branches in SIMMAP format 
# (see read.simmap()) and data for 1 or more continuous characters
# Written by Liam J. Revell 2010, 2011, 2013

evol.vcv<-function(scm.tre,dat,maxit=2000,vars=FALSE){

	## internal functions

	# function converts symmetric matrix to vector by taking upper diagonal
	to.upper<-function(mat){
		v<-vector(mode="numeric"); k<-1
		for(i in 1:nrow(mat)) for(j in i:ncol(mat)){
			v[k]<-mat[i,j]
			k<-k+1
		}
		return(v)
	}

	# function converts vector to upper diagonal matrix
	upper.diag<-function(vect){
		m<-(-1+sqrt(1+8*length(vect)))/2
		mat<-matrix(0,m,m); k<-1
		for(i in 1:m) for(j in i:m){
			mat[i,j]<-vect[k]
			k<-k+1
		}
		return(mat)
	}

	# function converts vector to symmetric matrix
	to.symmetric<-function(vect){
		m<-(-1+sqrt(1+8*length(vect)))/2
		mat<-matrix(0,m,m); k<-1
		for(i in 1:m) for(j in i:m){
			mat[i,j]<-vect[k]
			mat[j,i]<-mat[i,j]
			k<-k+1
		}
		return(mat)
	}

	# log-likelihood single rate matrix
	lik<-function(R.vect,y,D,a,C,n,m1){
		R<-to.symmetric(R.vect)
		logL1<--t(y-D%*%t(a))%*%solve(kronecker(R,C))%*%(y-D%*%t(a))/2-n*m*log(2*pi)/2-determinant(kronecker(R,C))$modulus[1]/2
	}

	# compute the log-likelihood (from the cholesky matrices)
	likelihood.cholR<-function(theta,y,C,D){
		m<-length(y)/dim(C)[1]; n<-length(y)/m; p<-dim(C)[3]
		cholR<-array(data=0,dim=c(m,m,p)); l<-1
		for(i in 1:p) for(j in 1:m) for(k in j:m){ 
			cholR[j,k,i]<-theta[l]; l<-l+1
		}
		V<-matrix(0,nrow(D),nrow(D))
		for(i in 1:p) V<-V+kronecker(t(cholR[,,i])%*%cholR[,,i],C[,,i])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# compute the log-likelihood (from the original matrices)
	likelihood.R<-function(theta,y,C,D){
		m<-length(y)/dim(C)[1]; n<-length(y)/m; p<-dim(C)[3]
		R<-array(data=0,dim=c(m,m,p)); l<-1
		for(i in 1:p) for(j in 1:m) for(k in j:m){ 
			R[j,k,i]<-theta[l]; R[k,j,i]<-theta[l]
			l<-l+1
		}
		V<-matrix(0,nrow(D),nrow(D))
		for(i in 1:p) V<-V+kronecker(R[,,i],C[,,i])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	## end internal functions

	# bookkeeping
	X<-as.matrix(dat)
	n<-nrow(X) # number of species
	m<-ncol(X) # number of traits
	p<-ncol(scm.tre$mapped.edge) # number of states
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	one<-matrix(1,n,1)
	y<-as.matrix(as.vector(X))

	# first compute C for the whole scm.tre
	C<-vcv.phylo(scm.tre)
	C<-C[rownames(dat),rownames(dat)]

	# find the ancestral state vector
	a<-colSums(solve(C))%*%X/sum(solve(C))

	# compute the MLE of R
	R<-t(X-one%*%a)%*%solve(C)%*%(X-one%*%a)/n

	# compute the log-likelihood
	logL1<-lik(to.upper(R),y,D,a,C,n,m)
	
	if(vars){
		# compute Hessian
		H<-hessian(lik,to.upper(R),y=y,D=D,a=a,C=C,n=n,m1=m)
		# compute variances
		vars.single<-to.symmetric(diag(solve(-H)))
		dimnames(vars.single)<-list(colnames(X),colnames(X))
	}

	# compute separate C for each state
	multi.tre<-list(); class(multi.tre)<-"multiPhylo"
	C<-array(dim=c(nrow(C),ncol(C),ncol(scm.tre$mapped.edge)))
	for(i in 1:ncol(scm.tre$mapped.edge)){
		multi.tre[[i]]<-scm.tre
		multi.tre[[i]]$edge.length<-scm.tre$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(scm.tre$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		C[,,i]<-temp[rownames(dat),rownames(dat)]
	}

	# compute the starting parameter values
	starting<-vector(mode="numeric")
	for(i in 1:p)
		starting<-c(starting,to.upper(chol(R)))

	# optimize using generic optimizer
	r=optim(starting,fn=likelihood.cholR,y=y,C=C,D=D,control=list(maxit=maxit))

	# convert parameter estimates to matrices
	R.i<-array(dim=c(m,m,p)); states<-vector(mode="character")
	for(i in 1:p){
		R.i[,,i]<-matrix(data=t(upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]))%*%upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]),m,m,dimnames=list(colnames(dat),colnames(dat)))
		states[i]=multi.tre[[i]]$state
	}
	dimnames(R.i)<-list(rownames(R),colnames(R),states)

	# log-likelihood for the multi-matrix model
	logL2<--r$value

	# if computing variances
	if(vars){
		# convert R.i to a vector
		R.vect<-vector(mode="numeric")
		for(i in 1:p)
			R.vect<-c(R.vect,to.upper(R.i[,,i]))
		# evaluate the Hessian
		H<-hessian(likelihood.R,R.vect,y=y,C=C,D=D)
		# convert diagonal of solve(H) to set of matrices
		Vars<-array(dim=c(m,m,p))
		Vh<-diag(solve(H))
		for(i in 1:p){
			Vars[,,i]<-matrix(data=t(to.symmetric(Vh[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))])),m,m)
			states[i]=multi.tre[[i]]$state
		}
		dimnames(Vars)<-list(rownames(R),colnames(R),states)
	}

	print(r)
	print(Vars)	
		
	# check for "false convergence"
	while(logL1>logL2||any(diag(Vh)<0)){
		message("Optimization may not have converged.  Changing starting value and repeating.")
		r=optim(starting+rnorm(n=length(starting),sd=0.1*mean(starting)),fn=likelihood.cholR,y=y,C=C,D=D,control=list(maxit=maxit))
		R.i<-array(dim=c(m,m,p)); states<-vector(mode="character")
		for(i in 1:p){
			R.i[,,i]<-matrix(data=t(upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]))%*%upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]),m,m,dimnames=list(colnames(dat),colnames(dat)))
			states[i]=multi.tre[[i]]$state
		}
		dimnames(R.i)<-list(rownames(R),colnames(R),states)
		logL2<--r$value
		if(vars){
			R.vect<-vector(mode="numeric")
			for(i in 1:p)
				R.vect<-c(R.vect,to.upper(R.i[,,i]))
			H<-hessian(likelihood.R,R.vect,y=y,C=C,D=D)
			Vars<-array(dim=c(m,m,p))
			Vh<-diag(solve(H))
			for(i in 1:p){
				Vars[,,i]<-matrix(data=t(to.symmetric(Vh[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))])),m,m)
				states[i]=multi.tre[[i]]$state
			}
			dimnames(Vars)<-list(rownames(R),colnames(R),states)
		}
		print(r)
		print(Vars)			
	}

	# report convergence
	if(r$convergence==0)
		converged<-"Optimization has converged."
	else
		converged<-"Optimization may not have converged.  Consider increasing maxit."

	# return results
	if(vars)
		return(list(R.single=R,vars.single=vars.single,logL1=as.numeric(logL1),k1=m*(m+1)/2+m,R.multiple=R.i,vars.multiple=Vars,logL.multiple=logL2,k2=p*m*(m+1)/2+m,P.chisq=pchisq(2*(logL2-as.numeric(logL1)),(p-1)*m*(m+1)/2,lower.tail=FALSE),convergence=converged))
	else
		return(list(R.single=R,logL1=as.numeric(logL1),k1=m*(m+1)/2+m,R.multiple=R.i,logL.multiple=logL2,k2=p*m*(m+1)/2+m,P.chisq=pchisq(2*(logL2-as.numeric(logL1)),(p-1)*m*(m+1)/2,lower.tail=FALSE),convergence=converged))

}
