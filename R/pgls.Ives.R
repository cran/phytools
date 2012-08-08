# implements method of Ives et al. 2007 for PGLS regression with sampling error
# written by Liam J. Revell 2012

pgls.Ives<-function(tree,X,y,Vx,Vy,Cxy){
	
	# likelihood function
	lik<-function(theta,C,x,y,Mx,My,Mxy){
		sig2x<-theta[1]; sig2y<-theta[2]; b1<-theta[3]
		a<-theta[4:5]
		n<-nrow(C)
		Psi<-matrix(0,2*n,2*n)
		Psi[1:n,1:n]<-sig2x*C+diag(Mx)
		Psi[n+1:n,1:n]<-Psi[1:n,n+1:n]<-b1*sig2x*C+diag(Mxy)
		Psi[n+1:n,n+1:n]<-b1^2*sig2x*C+sig2y*C+diag(My)
		z<-c(X,y)
		D<-kronecker(diag(rep(1,2)),matrix(rep(1,n)))
		L<--2*n/2*log(2*pi)-(1/2)*determinant(Psi,logarithm=TRUE)$modulus[1]-(1/2)*t(z-D%*%a)%*%solve(Psi)%*%(z-D%*%a)
		return(-L)
	}

	# perform calculation & organization
	C<-vcv.phylo(tree)
	X<-X[tree$tip.label]
	y<-y[tree$tip.label]
	Vx<-Vx[tree$tip.label]
	Vy<-Vy[tree$tip.label]
	Cxy<-Cxy[tree$tip.label]
	
	# get some reasonable starting values for optimization
	b<-runif(n=1,min=0,max=2)*lm(pic(X,tree)~pic(y,tree))$coefficients[2]; names(b)<-NULL
	sig2x<-runif(n=1,min=0,max=2)*mean(pic(X,tree)^2)
	sig2y<-runif(n=1,min=0,max=2)*mean(pic(y,tree)^2)
	a<-runif(n=2,min=-1,max=1)*c(mean(X),mean(y))

	# optimize regression model
	r<-optim(c(sig2x,sig2y,b,a),lik,C=C,x=X,y=y,Mx=Vx,My=Vy,Mxy=Cxy,method="L-BFGS-B",lower=c(0.0001,0.0001,-Inf,-Inf,-Inf),control=list(factr=1e10))

	# return r
	return(list(beta=c(r$par[5]-r$par[3]*r$par[4],r$par[3]),sig2x=r$par[1],sig2y=r$par[2],a=r$par[4:5],logL=-r$value))
}
