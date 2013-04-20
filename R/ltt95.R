# 95% CI on ltts
# written by Liam J. Revell 2013

ltt95<-function(trees,alpha=0.05,log=FALSE,method=c("lineages","times"),mode=c("median","mean"),...){
	method<-method[1]
	mode=mode[1]
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	X<-ltt(trees,plot=FALSE,gamma=FALSE)
	if(method=="times"){
		N<-length(X)
		tt<-sapply(X,function(x) max(x$times))
		zz<-max(tt)-tt
		for(i in 1:N) X[[i]]$times<-X[[i]]$times+zz[i]
		n<-sapply(X,function(x) max(x$ltt))
		if(all(n==max(n))) n<-max(n) 
		else stop("for method='times' all trees must contain the same numer of lineages")
		LL<-sapply(X,function(x) x$times[1:length(x$times)])
		ii<-floor(alpha/2*N)
		jj<-ceiling((1-alpha/2)*N)
		low<-apply(LL,1,function(x) sort(x)[ii])
		high<-apply(LL,1,function(x) sort(x)[jj])
		ll<-if(mode=="median") apply(LL,1,function(x) median(x)[1]) else rowMeans(LL)
		plot(ll,c(1:n,n),lwd=2,xlim=c(0,max(tt)),
			type=if(mode=="median") "s" else "l",main=NULL,
			xlab="time",ylab="lineages",log=if(log) "y" else "")
		lines(low,c(1:n,n),lty="dashed",type=if(mode=="median") "s" else "l")
		lines(high,c(1:n,n),lty="dashed",type=if(mode=="median") "s" else "l")
		XX<-cbind(c(1:n,n),low,ll,high)
		colnames(XX)<-c("lineages","low(time)","time","high(time)")
		rownames(XX)<-NULL
	} else if(method=="lineages"){
		N<-length(X)
		tt<-sapply(X,function(x) max(x$times))
		zz<-max(tt)-tt
		for(i in 1:N) X[[i]]$times<-X[[i]]$times+zz[i]
		tt<-0:res*max(tt)/res
		ll<-low<-high<-vector()
		for(i in 1:length(tt)){
			ss<-vector()
			for(j in 1:N){
				ii<-2
				while(tt[i]>X[[j]]$times[ii]&&ii<length(X[[j]]$times)) ii<-ii+1
				ss[j]<-X[[j]]$ltt[ii-1]
			}
			ll[i]<-if(mode=="median") median(ss) else mean(ss)
			low[i]<-sort(ss)[floor(alpha/2*N)]
			high[i]<-sort(ss)[ceiling((1-alpha/2)*N)]
		}
		plot(tt,ll,ylim=c(min(low),max(high)),lwd=2,
			type=if(mode=="median") "s" else "l",main=NULL,
			xlab="time",ylab="lineages",log=if(log) "y" else "")
		lines(tt,low,lty="dashed",type="s")
		lines(tt,high,lty="dashed",type="s")
		XX<-cbind(tt,low,ll,high)
		colnames(XX)<-c("time","low(lineages)","lineages","high(lineages)")
		rownames(XX)<-NULL
	}
	invisible(XX)
}