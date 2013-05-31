# function to plot probability or trait value by branch
# written by Liam J. Revell 2013

plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
	mode<-mode[1]
	if(mode=="tips"){
		x<-c(x[tree$tip.label],fastAnc(tree,x))
		names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
		XX<-matrix(x[tree$edge],nrow(tree$edge),2)
		x<-rowMeans(XX)
	} else if(mode=="nodes"){
		XX<-matrix(x[tree$edge],nrow(tree$edge),2)
		x<-rowMeans(XX)
	}
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-6
	if(palette=="heat.colors") cols<-heat.colors(n=1000)
	if(palette=="gray") cols<-gray(1000:1/1000)
	if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
	if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
	breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
	whichColor<-function(p,cols,breaks){
		i<-1
		while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
		cols[i]
	}
	colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
	par(lend=2)
	plot.phylo(tree,no.margin=TRUE,edge.width=4,edge.color=colors,label.offset=0.01*max(nodeHeights(tree)),lend=2,new=FALSE,...)
	if(legend==TRUE&&is.logical(legend)) legend<-round(0.3*max(nodeHeights(tree)),2)
	if(legend){
		if(hasArg(title)) title<-list(...)$title
		else title<-NULL
		if(hasArg(digits)) digits<-list(...)$digits
		else digits<-1
		add.color.bar(legend,cols,title,xlims,digits,prompt=FALSE)
	}
}

# function to add color bar
# written by Liam J. Revell 2013

add.color.bar<-function(leg,cols,title=NULL,lims=c(0,1),digits=1,prompt=TRUE,...){
	if(prompt){
		cat("Click where you want to draw the bar\n")
		x<-unlist(locator(1))
		y<-x[2]
		x<-x[1]
	} else {
		if(hasArg(x)) x<-list(...)$x
		else x<-0
		if(hasArg(y)) y<-list(...)$y
		else y<-0
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	X<-x+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
	Y<-cbind(rep(y,length(cols)),rep(y,length(cols))) 		
	lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=4+2,lend=2) 
	for(i in 1:length(cols)) lines(X[i,],Y[i,],col=cols[i],lwd=4,lend=2)
	text(x=x,y=y,round(lims[1],digits),pos=3,cex=fsize)
	text(x=x+leg,y=y,round(lims[2],digits),pos=3,cex=fsize)
	if(is.null(title)) title<-"P(state=1)"
	text(x=(2*x+leg)/2,y=y,title,pos=3,cex=fsize)
	text(x=(2*x+leg)/2,y=y,paste("length=",round(leg,3),sep=""),pos=1,cex=fsize)
}
