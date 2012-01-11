# function is simplified version of evol.vcv
# written by Liam J. Revell 2011

evolvcv.lite<-function(tree,X,maxit=2000,tol=1e-10){

	# model 1: common variances & correlation
	lik1<-function(theta,C,D,y){
		v<-theta[1:2]; r<-theta[3]
		R<-matrix(c(v[1],r*sqrt(v[1]*v[2]),r*sqrt(v[1]*v[2]),v[2]),2,2)
		V<-kronecker(R,C)
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 2: different variances, same correlation
	lik2<-function(theta,C,D,y){
		v1<-theta[1:2]; v2<-theta[3:4]; r<-theta[5]
		R1<-matrix(c(v1[1],r*sqrt(v1[1]*v1[2]),r*sqrt(v1[1]*v1[2]),v1[2]),2,2)
		R2<-matrix(c(v2[1],r*sqrt(v2[1]*v2[2]),r*sqrt(v2[1]*v2[2]),v2[2]),2,2)
		V<-kronecker(R1,C[[1]])+kronecker(R2,C[[2]])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 3: same variances, different correlation
	lik3<-function(theta,C,D,y){
		v<-theta[1:2]; r1<-theta[3]; r2<-theta[4]
		R1<-matrix(c(v[1],r1*sqrt(v[1]*v[2]),r1*sqrt(v[1]*v[2]),v[2]),2,2)
		R2<-matrix(c(v[1],r2*sqrt(v[1]*v[2]),r2*sqrt(v[1]*v[2]),v[2]),2,2)
		V<-kronecker(R1,C[[1]])+kronecker(R2,C[[2]])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 4: everything different
	lik4<-function(theta,C,D,y){
		v1<-theta[1:2]; v2<-theta[3:4]; r1<-theta[5]; r2<-theta[6]
		R1<-matrix(c(v1[1],r1*sqrt(v1[1]*v1[2]),r1*sqrt(v1[1]*v1[2]),v1[2]),2,2)
		R2<-matrix(c(v2[1],r2*sqrt(v2[1]*v2[2]),r2*sqrt(v2[1]*v2[2]),v2[2]),2,2)
		V<-kronecker(R1,C[[1]])+kronecker(R2,C[[2]])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}
		
	# done internal functions

	# bookkeeping
	X<-as.matrix(X)
	X<-X[tree$tip.label,]
	n<-nrow(X) # number of species
	m<-ncol(X) # number of traits
	if(m!=2) stop("number of traits must equal 2")
	p<-ncol(tree$mapped.edge) # number of states
	if(p!=2) stop("number of states must equal 2")
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	y<-as.matrix(as.vector(X))
	mC<-multiC(tree)
	C<-mC[[1]]+mC[[2]]
	sv1<-mean(pic(X[,1],tree)^2)
	sv2<-mean(pic(X[,2],tree)^2)
	sr<-mean(pic(X[,1],tree)*pic(X[,2],tree))/sqrt(sv1*sv2)

	# now optimize models
	res1<-optim(runif(3)*c(sv1,sv2,sr),lik1,C=C,D=D,y=y,method="L-BFGS-B",lower=tol+c(0,0,-1),upper=c(Inf,Inf,1)-tol)
	res2<-optim(runif(5)*c(sv1,sv2,sv1,sv2,sr),lik2,C=mC,D=D,y=y,method="L-BFGS-B",lower=tol+c(0,0,0,0,-1),upper=c(Inf,Inf,Inf,Inf,1)-tol)
	res3<-optim(runif(4)*c(sv1,sv2,sr,sr),lik3,C=mC,D=D,y=y,method="L-BFGS-B",lower=tol+c(0,0,-1,-1),upper=c(Inf,Inf,1,1)-tol)
	res4<-optim(runif(6)*c(sv1,sv2,sv1,sv2,sr,sr),lik4,C=mC,D=D,y=y,method="L-BFGS-B",lower=tol+c(0,0,0,0,-1,-1),upper=c(Inf,Inf,Inf,Inf,1,1)-tol)

	m1<-list(description="common rates, common correlation",
				R=matrix(c(res1$par[1],res1$par[3]*sqrt(res1$par[1]*res1$par[2]),res1$par[3]*sqrt(res1$par[1]*res1$par[2]),res1$par[2]),2,2),
				logLik=-res1$value,
				convergence=res1$convergence,
				k=length(res1$par)+1,
				AIC=2*(length(res1$par)+1)+2*res1$value)
	m2<-list(description="different rates, common correlation",
				R1=matrix(c(res2$par[1],res2$par[5]*sqrt(res2$par[1]*res2$par[2]),res2$par[5]*sqrt(res2$par[1]*res2$par[2]),res2$par[2]),2,2),
				R2=matrix(c(res2$par[3],res2$par[5]*sqrt(res2$par[3]*res2$par[4]),res2$par[5]*sqrt(res2$par[3]*res2$par[4]),res2$par[4]),2,2),
				logLik=-res2$value,
				convergence=res2$convergence,
				k=length(res2$par)+1,
				AIC=2*(length(res2$par)+1)+2*res2$value)
	m3<-list(description="common rates, different correlation",
				R1=matrix(c(res3$par[1],res3$par[3]*sqrt(res3$par[1]*res3$par[2]),res3$par[3]*sqrt(res3$par[1]*res3$par[2]),res3$par[2]),2,2),
				R2=matrix(c(res3$par[1],res3$par[4]*sqrt(res3$par[1]*res3$par[2]),res3$par[4]*sqrt(res3$par[1]*res3$par[2]),res3$par[2]),2,2),
				logLik=-res3$value,
				convergence=res3$convergence,
				k=length(res3$par)+1,
				AIC=2*(length(res3$par)+1)+2*res3$value)
	m4<-list(description="no common structure",
				R1=matrix(c(res4$par[1],res4$par[5]*sqrt(res4$par[1]*res4$par[2]),res4$par[5]*sqrt(res4$par[1]*res4$par[2]),res4$par[2]),2,2),
				R2=matrix(c(res4$par[3],res4$par[6]*sqrt(res4$par[3]*res4$par[4]),res4$par[6]*sqrt(res4$par[3]*res4$par[4]),res4$par[4]),2,2),
				logLik=-res4$value,
				convergence=res4$convergence,
				k=length(res4$par)+1,
				AIC=2*(length(res4$par)+1)+2*res4$value)
	
	return(list(model1=m1,model2=m2,model3=m3,model4=m4))

}

# function
# written by Liam J. Revell

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
