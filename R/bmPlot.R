# visualize discrete time Brownian simulation on a tree
# written by Liam J. Revell 2012

bmPlot<-function(tree,type="BM",anc=0,sig2=1/1000,ngen=1000,...){
	tr<-reorder(tree)
	H<-nodeHeights(tr)
	tr$edge.length<-round(tr$edge.length/max(H)*ngen)
	h<-nodeHeights(tr)[,1]

	bmSim<-function(start,n)
		cumsum(c(start,rnorm(n,sd=sqrt(sig2))))

	X<-T<-list()
	N<-length(tr$tip)
	for(i in 1:nrow(tr$edge)){
		if(tr$edge[i,1]==(N+1)) X[[i]]<-bmSim(anc,tr$edge.length[i])
		else {
			parent<-match(tr$edge[i,1],tr$edge[,2])
			X[[i]]<-bmSim(X[[parent]][length(X[[parent]])],tr$edge.length[i])
		}
		T[[i]]<-h[i]+0:tr$edge.length[i]
	}
	
	if(type=="BM") bm(X,T)
	else if(type=="threshold") th(X,T,...)

	Y<-matrix(NA,nrow(tr$edge),2)
	Y[,1]<-sapply(X,function(x) x[1])
	Y[,2]<-sapply(X,function(x) x[length(x)])

	x<-Y[tr$edge[,2]%in%1:N,2]
	names(x)<-tr$tip.label[tr$edge[tr$edge[,2]%in%1:N,2]]
	a<-c(anc,Y[tr$edge[,2]%in%(N+2:tr$Nnode),2])
	names(a)<-c(N+1,tr$edge[tr$edge[,2]%in%(N+2:tr$Nnode),2])
	
	return(c(x[tree$tip],a[as.character(N+1:tree$Nnode)]))
}

# plots type="BM"
bm<-function(X,T){
	minX<-min(sapply(X,min))
	maxX<-max(sapply(X,max))
	plot(X[[1]][1],T[[1]][1],ylim=c(0,max(sapply(T,max))),xlim=c(minX,maxX),ylab="time",xlab="phenotype")
	for(i in 1:length(X)) lines(X[[i]],T[[i]])
}

# plots type="threshold"
th<-function(X,T,...){
	minX<-min(sapply(X,min))
	maxX<-max(sapply(X,max))
	if(hasArg(thresholds)) thresholds<-list(...)$thresholds
	else stop("no thresholds provided for type='threshold'")
	if(hasArg(cols)) cols<-list(...)$cols
	else {
		cols<-c("black","red","blue","green")
		while(length(thresholds)>(length(cols)-1)) cols<-c(cols,c("black","red","blue","green"))
	}
	thresholds<-thresholds+1e-16
	plot(X[[1]][1],T[[1]][1],ylim=c(0,max(sapply(T,max))),xlim=c(minX,maxX),ylab="time",xlab="liability")
	for(i in 1:length(X))
		if(length(X[[i]])>1)
			for(j in 1:(length(X[[i]])-1))
				if(threshCol(X[[i]][j],thresholds,cols)==threshCol(X[[i]][j+1],thresholds,cols))
					lines(X[[i]][c(j,j+1)],T[[i]][c(j,j+1)],col=threshCol(X[[i]][j+1],thresholds,cols))
				else {
					t<-thresholds[which(sort(c(thresholds,min(X[[i]][c(j,j+1)])))==min(X[[i]][c(j,j+1)]))]
					lines(c(X[[i]][j],t),c(T[[i]][j],T[[i]][j]+abs(diff(c(X[[i]][j],t))/diff(X[[i]][c(j,j+1)]))),col=threshCol(X[[i]][j],thresholds,cols))
					lines(c(X[[i]][j+1],t),c(T[[i]][j+1],T[[i]][j+1]-abs(diff(c(X[[i]][j+1],t))/diff(X[[i]][c(j,j+1)]))),col=threshCol(X[[i]][j+1],thresholds,cols))
				}
	for(i in 1:length(thresholds)) lines(rep(thresholds[i],2),c(0,max(sapply(T,max))),lty="dashed")
}

# returns a color based on position relative to thresholds
threshCol<-function(x,thresholds,cols){
	t<-c(-Inf,thresholds,Inf); i<-1
	while(x>t[i]) i<-i+1
	return(cols[i-1])
}

