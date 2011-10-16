# simulates with multiple evolutionary rates in different parts of the tree
# written by Liam J. Revell 2011

sim.rates<-function(mtree,sig2,anc=0,nsim=1,plot=F){
	if(class(mtree)!="phylo") stop("mtree should be an object of class 'phylo'")
	if(is.null(mtree$mapped.edge)){
		message("mtree does not contain a mapped discrete character history, using fastBM")
		X<-fastBM(mtree,a=anc,sig2=sig2[1],nsim=nsim)
	} else {
		# first name (if necessary) and reorder sig2
		if(is.null(names(sig2))){
			message("names absent from sig2: assuming same order as $mapped.edge")
			if(length(sig2)==ncol(mtree$mapped.edge)) names(sig2)<-colnames(mtree$mapped.edge)
			else stop("the number of elements in sig2 should match the number of rows in mapped.edge")
		}
		sig2<-sig2[colnames(mtree$mapped.edge)]
		# now create a tree for simulation
		edge.length<-rep(0,nrow(mtree$edge))
		# scale the edge lengths by the rate
		for(i in 1:ncol(mtree$mapped.edge))
			edge.length<-edge.length+sig2[i]*mtree$mapped.edge[,i]
		names(edge.length)<-NULL
		tree<-list(Nnode=mtree$Nnode,edge=mtree$edge,tip.label=mtree$tip.label,edge.length=edge.length)
		class(tree)<-"phylo"
		if(plot) plot(tree)
		# simulate
		X<-fastBM(tree,a=anc,nsim=nsim)
	}
	return(X)
}
