# this function computes the data for a lineage-through-time plot and (optionally) creates this plot
# the function does not require a tree that is ultrametric.  Optionally, the function can remove extinct
# species from the phylogeny
# if the input tree is an object of class "multiPhylo" then the function will simultaneously plot all ltts
# written by Liam J. Revell 2010, 2011, 2012

ltt<-function(tree,plot=TRUE,drop.extinct=FALSE,log.lineages=TRUE,gamma=TRUE){
	# set tolerance
	tol<-1e-6
	# check 'phylo' object
	if(class(tree)!="phylo"&&class(tree)!="multiPhylo") stop("tree object must be of class 'phylo.'")
	if(class(tree)=="multiPhylo"){
		X<-lapply(tree,ltt,plot=FALSE,drop.extinct=drop.extinct,log.lineages=log.lineages,gamma=gamma)
		max.lineages<-max(sapply(X,function(y) max(y$ltt)))
		max.time<-max(sapply(X,function(y) max(y$times)))
		if(plot==TRUE){
			if(log.lineages==TRUE) plot(X[[1]]$time,log(X[[1]]$ltt),"s",xlab="time",ylab="log(lineages)",xlim=c(0,max.time),ylim=c(0,log(max.lineages)))
			else plot(X[[1]]$time,X[[1]]$ltt,"s",xlab="time",ylab="lineages",xlim=c(0,max.time),ylim=c(0,max.lineages))
			for(i in 2:length(trees))
				if(log.lineages==TRUE) lines(X[[i]]$time,log(X[[i]]$ltt),"s",xlab="time",ylab="log(lineages)",xlim=c(0,max.time),ylim=c(0,log(max.lineages)))
				else lines(X[[i]]$time,X[[i]]$ltt,"s",xlab="time",ylab="lineages",xlim=c(0,max.time),ylim=c(0,max.lineages))
		}
		return(X)
	} else {
		# reorder the tree
		tree<-reorder.phylo(tree,order="cladewise")
		# internal functions
		# drop extinct tips
		drop.extinct.tips<-function(phy){
			temp<-diag(vcv(phy))
			if(length(temp[temp<(max(temp)-tol)])>0)
				pruned.phy<-drop.tip(phy,names(temp[temp<(max(temp)-tol)]))
			else
				pruned.phy<-phy
			return(pruned.phy)
		}
		# first, if drop.extinct==TRUE
		if(drop.extinct==TRUE)
			tree<-drop.extinct.tips(tree)
		# compute node heights
		root<-length(tree$tip)+1
		node.height<-matrix(NA,nrow(tree$edge),2)
		for(i in 1:nrow(tree$edge)){
			if(tree$edge[i,1]==root){
				node.height[i,1]<-0.0
				node.height[i,2]<-tree$edge.length[i]
			} else {
				node.height[i,1]<-node.height[match(tree$edge[i,1],tree$edge[,2]),2]
				node.height[i,2]<-node.height[i,1]+tree$edge.length[i]
			}
		}
		ltt<-vector()
		tree.length<-max(node.height) # tree length
		n.extinct<-sum(node.height[tree$edge[,2]<=length(tree$tip),2]<(tree.length-tol))
		# fudge things a little bit
		node.height[tree$edge[,2]<=length(tree$tip),2]<-node.height[tree$edge[,2]<=length(tree$tip),2]+1.1*tol
		time<-c(0,node.height[,2]); names(time)<-as.character(c(root,tree$edge[,2]))
		temp<-vector()
		time<-time[order(time)]
		time<-time[1:(tree$Nnode+n.extinct+1)] # times
		# get numbers of lineages
		for(i in 1:(length(time)-1)){
			ltt[i]<-0
			for(j in 1:nrow(node.height))
				ltt[i]<-ltt[i]+(time[i]>=(node.height[j,1]-tol)&&time[i]<=(node.height[j,2]-tol))
		}
		ltt[i+1]<-0
		for(j in 1:nrow(node.height))
			ltt[i+1]<-ltt[i+1]+(time[i+1]<=(node.height[j,2]+tol))
		# set names (these are the node indices)
		names(ltt)<-names(time)
		# append 0,1 for time 0
		ltt<-c(1,ltt); time<-c(0,time)
		# subtract fudge factor
		time[length(time)]<-time[length(time)]-1.1*tol
		# plot ltt
		if(plot==TRUE){
			if(log.lineages==TRUE)
				plot(time,log(ltt),"s",xlab="time",ylab="log(lineages)")
			else
				plot(time,ltt,"s",xlab="time",ylab="lineages")
		}
		if(gamma==FALSE)
			return(list(ltt=ltt,times=time))
		else{
			gam<-gammatest(list(ltt=ltt,times=time))
			return(list(ltt=ltt,times=time,gamma=gam$gamma,p=gam$p))
		}
	}
}

# function computes the gamma-statistic & a two-tailed P-value
# written by Liam Revell 2011

gammatest<-function(x){
	n<-max(x$ltt)
	g<-vector()
	for(i in 2:(length(x$times))) g[i-1]<-x$times[i]-x$times[i-1]
	T<-sum((2:n)*g[2:n])
	doublesum<-0
	for(i in 1:(n-1)) for(k in 1:i) doublesum<-doublesum+k*g[k]
	gamma<-(1/(n-2)*doublesum-T/2)/(T*sqrt(1/(12*(n-2))))
	p<-2*pnorm(abs(gamma),lower.tail=F)
	return(list(gamma=gamma,p=p))
}
