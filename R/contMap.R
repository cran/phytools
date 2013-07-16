# function plots reconstructed values for ancestral characters along the edges of the tree
# written by Liam J. Revell 2012-2013

contMap<-function(tree,x,res=100,fsize=NULL,ftype=NULL,lwd=4,legend=NULL,lims=NULL,outline=TRUE,sig=3,type="phylogram",plot=TRUE,...){
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	h<-max(nodeHeights(tree))
	steps<-0:res/res*max(h)
	H<-nodeHeights(tree)
	a<-fastAnc(tree,x)
	y<-c(a,x[tree$tip.label]); names(y)[1:length(tree$tip)+tree$Nnode]<-1:length(tree$tip)
	A<-matrix(y[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
	cols<-rainbow(1001,start=0,end=0.7); names(cols)<-0:1000
	if(is.null(lims)) lims<-c(min(x),max(x))
	trans<-0:1000/1000*(lims[2]-lims[1])+lims[1]; names(trans)<-0:1000
	for(i in 1:nrow(tree$edge)){
		XX<-cbind(c(H[i,1],steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))]),
			c(steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))],H[i,2]))-H[i,1]
		YY<-rowMeans(XX)
		b<-vector()
		for(j in 1:length(YY))
			b[j]<-(A[i,1]/YY[j]+A[i,2]/(max(XX)-YY[j]))/(1/YY[j]+1/(max(XX)-YY[j]))
		d<-sapply(b,getState,trans=trans)
		tree$maps[[i]]<-XX[,2]-XX[,1]
		names(tree$maps[[i]])<-d
	}
	xx<-list(tree=tree,cols=cols,lims=lims)
	class(xx)<-"contMap"
	if(plot) plot.contMap(xx,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,sig=sig,type=type,mar=mar)
	invisible(xx)
}

# function
# written by Liam J. Revell 2012

getState<-function(x,trans){
	i<-1
	while(x>trans[i]){
		state<-names(trans)[i]
		i<-i+1
	}
	return(state)
}

# function
# written by Liam J. Revell 2012-2013

plot.contMap<-function(x,...){
	if(class(x)=="contMap"){
		lims<-x$lims
		x<-list(tree=x$tree,cols=x$cols)
		class(x)<-"densityMap"
	} else stop("x should be an object of class 'contMap'")
	H<-nodeHeights(x$tree)
	# get & set optional arguments
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-TRUE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-4
	if(hasArg(sig)) sig<-list(...)$sig
	else sig<-3
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	if(is.null(legend)) legend<-0.5*max(H)
	if(is.null(fsize)) fsize<-c(1,1)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(is.null(ftype)) ftype<-c("i","reg")
	if(length(ftype)==1) ftype<-c(ftype,"reg")
	# done optional arguments
	leg.txt<-c(round(lims[1],sig),"trait value",round(lims[2],sig))
	plot(x,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,leg.txt=leg.txt,type=type,mar=mar)
}

