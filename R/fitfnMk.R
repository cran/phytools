fitfnMk<-function(tree,x,model="polynomial",degree=2,...){
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(model!="polynomial"){
		stop("Sorry. Only available model (so far) is \"polynomial\". Stopping.\n")
	}
	if(hasArg(niter)) niter<-list(...)$niter
	else niter<-10
	if(niter>1){
		if(hasArg(parallel)) parallel<-list(...)$parallel
		else parallel<-FALSE
		if(hasArg(opt.method)) opt.method<-list(...)$opt.method
		else opt.method<-"nlminb"
	} else parallel<-FALSE
	if(niter>1){
		args<-list(...)
		args$niter<-1
		args$tree<-tree
		args$x<-x
		args$degree<-degree
		if(parallel&&(opt.method!="DEoptim")){
			if(hasArg(ncores)) ncores<-list(...)$ncores
			else ncores<-min(c(detectCores()-1,niter))
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
			if(!quiet){
				cat(paste("Opened cluster with",ncores,"cores.\n"))
				cat("Running optimization iterations in parallel.\n")
				cat("Please wait....\n")
				flush.console()
			}
			fits<-foreach(i = 1:niter)%dopar%{
				do.call(fitfnMk,args)
			}
			stopCluster(cl=mc)
			lnL<-sapply(fits,logLik)
		} else {
			fits<-list()
			for(i in 1:niter){ 
				fits[[i]]<-do.call(fitfnMk,args)
				lnL<-sapply(fits,logLik)
				if(!quiet){
					cat(paste("log-likelihood from current iteration:", 
						round(logLik(fits[[i]]),4),"\n"))
					cat(paste(" --- Best log-likelihood so far:", 
						round(max(lnL), 4),"---\n"))
						flush.console()
				}
			}
		}
		object<-fits[[which(lnL==max(lnL))[1]]]
		object$all.fits<-fits
	} else {
		if(hasArg(trace)) trace<-list(...)$trace
		else trace<-0
		if(hasArg(start)) start<-list(...)$start
		else start<-NULL
		if(hasArg(opt.method)) opt.method<-list(...)$opt.method
		else opt.method<-"nlminb"
		if(length(degree)==1) degree<-rep(degree,2)
		if(is.matrix(x)){
			levs<-colnames(x)
		} else if(is.numeric(x)){
			levs<-min(x):max(x)
			x<-to.matrix(x,levs)
		} else if(is.factor(x)){
			if(suppressWarnings(all(!is.na(as.numeric(levels(x)))))){
				levs<-min(as.numeric(levels(x))):max(as.numeric(levels(x)))
				x<-to.matrix(x,levs)
			} else {
				levs<-sort(levels(x))
				x<-to.matrix(x,levs)
			}
		} else if(is.character(x)){
			if(suppressWarnings(all(!is.na(as.numeric(x))))){
				levs<-min(as.numeric(x)):max(as.numeric(x))
				x<-to.matrix(x,levs)
			} else {
				levs<-sort(unique(x))
				x<-to.matrix(x,levs)
			}
		}
		x<-x[tree$tip.label,]
		k<-ncol(x)	
		if(hasArg(pi)) pi<-list(...)$pi
		else pi<-"equal"
		if(is.numeric(pi)) root.prior<-"given"
		if(pi[1]=="equal"){ 
			pi<-setNames(rep(1/k,k),levs)
			root.prior<-"flat"
		} else if(pi[1]=="fitzjohn"){ 
			root.prior<-"nuisance"
		} else if(pi[1]=="mle") root.prior<-"it's MLE"	
		lik<-function(par,pw,Y,pi=pi,degree=degree){
			k<-ncol(Y)
			m<-1:(k-1)-0.5
			q1<-rep(0,length(m))
			for(i in 0:degree[1]) q1<-q1+par[i+1]*m^(degree[1]-i)
			q2<-rep(0,length(m))
			for(i in 0:degree[2]) q2<-q2+par[degree[1]+i+2]*m^(degree[2]-i)
			q1[q1<0]<-0
			q2[q2<0]<-0
			if(all(q1<0)||all(q2<0)){ 
				return(Inf)
			} else {
				MODEL<-matrix(0,k,k,dimnames=list(colnames(Y),colnames(Y)))
				MODEL[cbind(1:(k-1),2:k)]<-1:(k-1)
				MODEL[cbind(2:k,1:(k-1))]<-k:(2*k-2)
				return(-pruning(c(q1,q2),pw,Y,model=MODEL,pi=pi))
			}
		}
		pw<-reorder(tree,"postorder")
		xx<-0:(k-2)+0.5
		if(is.null(start)) q1_start<-q2_start<--1
		else if(start[1]=="smart"){
			cat("Numerically optimizing simple equal-rates ordered model\n")
			cat("  to get better random starting values....\n\n")
			MODEL<-matrix(0,k,k,dimnames=list(colnames(x),colnames(x)))
			MODEL[cbind(1:(k-1),2:k)]<-1
			MODEL[cbind(2:k,1:(k-1))]<-2
			RATES<-fitMk(pw,x,model=MODEL,pi=pi,lik.func="pruning")$rates
			start<-rep(0,sum(degree)+2)
			start[c(degree[1]+1,sum(degree)+2)]<-RATES
			start<-start+runif(n=sum(degree)+2,min=-0.0001*mean(RATES),
				max=0.0001*mean(RATES))
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
			q1_start[q1_start<0]<-0
			q2_start[q2_start<0]<-0
			if(all(q1_start==0)&&all(q2_start==0)) q1_start<-q2_start<--1	
		} else if(is.numeric(start)) {
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
			q1_start[q1_start<0]<-0
			q2_start[q2_start<0]<-0
			if(all(q1_start==0)&&all(q2_start==0)) q1_start<-q2_start<--1
		} else q1_start<-q2_start<--1
		while(any(q1_start<0)||any(q2_start<0)){
			start<-runif(n=sum(degree)+2)
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
		}
		if(opt.method=="nlminb"){
			fit<-nlminb(start,lik,pw=pw,Y=x,pi=pi,degree=degree,
				control=list(trace=trace))
		} else if(opt.method=="DEoptim"){
			if(hasArg(maxit)) maxit<-list(...)$maxit
			else maxit<-1000
			if(hasArg(prompt)) prompt<-list(...)$prompt
			else prompt<-FALSE
			if(hasArg(parallel)) parallel<-list(...)$parallel
			else parallel<-FALSE
			if(hasArg(trace)) trace<-list(...)$trace
			else trace<-if(maxit>100) 100 else 1
			parallelType<-if(parallel) 1 else 0
			if(hasArg(NP)) NP<-list(...)$NP
			else NP<-100*length(start)
			if(hasArg(F)) F<-list(...)$F
			else F<-0.8
			if(hasArg(CR)) CR<-list(...)$CR
			else CR<-0.9
			initialpop<-matrix(rep(start,NP)*runif(NP*length(start)),
				NP,length(start),byrow=TRUE)
			if(parallel){
				if(hasArg(ncores)) ncores<-list(...)$ncores
				else ncores<-detectCores()-1
				mc<-makeCluster(ncores,type="PSOCK")
				registerDoParallel(cl=mc)
			}
			if(trace>0) cat("\n")
			fit<-DEoptim(lik,
				lower=rep(-10*max(start),length(start)),
				upper=rep(10*max(start),length(start)),
				pw=pw,Y=x,pi=pi,degree=degree,
				control=list(itermax=maxit,
				initialpop=initialpop,
				cluster=if(parallel) mc else NULL,
				trace=trace,F=F,CR=CR,NP=NP))
			if(prompt){ 
				cat("\nMax iterations (itermax) reached. Continue? (yes/no):")
				continue<-readline()
			} else continue<-"no"
			while(continue=="yes"&&prompt){
				fit<-DEoptim(lik,
					lower=rep(-10*max(start),length(start)),
					upper=rep(10*max(start),length(start)),
					pw=pw,Y=x,pi=pi,degree=degree,
					control=list(itermax=maxit,
					initialpop=fit$member$pop,
					cluster=if(parallel) mc else NULL,
					trace=trace,F=F,CR=CR,NP=NP))
				cat("\nMax iterations (itermax) reached. Continue? (yes/no):")
				continue<-readline()
			}
			if(trace>0) cat("\n")
			stopCluster(cl=mc)
			fit<-list(
				par=fit$optim$bestmem,
				objective=fit$optim$bestval,
				iterations=fit$optim$iter,
				evaluations=fit$optim$nfeval)
		}
		q1_est<-rep(0,length(xx))
		for(i in 0:degree[1]) q1_est<-q1_est+fit$par[i+1]*xx^(degree[1]-i)
		q2_est<-rep(0,length(xx))
		for(i in 0:degree[2]) q2_est<-q2_est+fit$par[degree[1]+i+2]*xx^(degree[2]-i)
		q1_est[q1_est<0]<-0
		q2_est[q2_est<0]<-0
		index.matrix<-matrix(0,k,k,dimnames=list(colnames(x),colnames(x)))
		index.matrix[cbind(1:(k-1),2:k)]<-1:(k-1)
		index.matrix[cbind(2:k,1:(k-1))]<-k:(2*k-2)
		lik.f<-function(par) -lik(par,pw=pw,Y=x,pi=pi,degree=degree)
		object<-list(
			logLik=-fit$objective,
			rates=c(q1_est,q2_est),
			index.matrix=index.matrix,
			states=levs,
			pi=pi,
			method=opt.method,
			root.prior=root.prior,
			opt_results=fit[c("convergence","iterations","evaluations","message")],
			par=fit$par,
			degree=degree,
			data=x,
			tree=tree,
			lik=lik.f)
		class(object)<-c("fitfnMk","fitMk")
	}	
	object
}

plot.fitfnMk<-function(x,...){
	k<-length(x$states)
	q1<-x$rates[1:(k-1)]
	q2<-x$rates[k:(2*k-2)]
	xx<-0:(k-2)+0.5
	plot(xx,q1,type="b",col="blue",bty="n",las=1,
		axes=FALSE,xlab="",ylab="transition rate (q)",
		ylim=c(0,max(c(q1,q2))))
	lines(xx,q2,type="b",col="red")
	labs<-mapply(function(x,y) bquote(.(x) %<->% .(y)),
		x=x$states[1:(k-1)],y=x$states[2:k])
	axis(1,at=seq(0.5,k-1.5,by=1),labels=rep("",k-1))
	nulo<-mapply(mtext,text=labs,at=seq(0.5,k-1.5,by=1),
		MoreArgs=list(side=1,line=1,las=3,cex=0.7))
	axis(2,las=1,cex.axis=0.8)
	grid()
	legend("bottomleft",c("forward","backward"),
		col=c("blue","red"),
		lty="solid",pch=1,cex=0.8)
}

logLik.fitfnMk<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-length(object$par)
	lik
}

## print method for objects of class "fitMk"
print.fitfnMk<-function(x,digits=6,...){
	cat("Object of class \"fitfnMk\".\n\n")
	cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-x$states
	print(round(Q,digits))
	cat("\nFitted (or set) value of pi:\n")
	print(round(x$pi,digits))
	cat(paste("due to treating the root prior as (a) ",x$root.prior,".\n",
		sep=""))
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",
		sep=""))
	if(!is.null(x$opt_results$convergence)){
		if(x$opt_results$convergence==0) 
			cat("R thinks it has found the ML solution.\n\n")
		else cat("R thinks optimization may not have converged.\n\n")
	}
}
