# function plots reconstructed values for ancestral characters along the edges of the tree
# written by Liam J. Revell 2012

contMap<-function(tree,x,res=100,fsize=NULL,ftype=NULL,lwd=4,legend=NULL,lims=NULL,outline=TRUE,sig=3){
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
	if(is.null(legend)) legend<-0.5*max(H)
	if(is.null(fsize)) fsize<-c(1,1)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(is.null(ftype)) ftype<-c("i","reg")
	if(length(ftype)==1) ftype<-c(ftype,"reg")
	if(legend){
		if(legend>max(H)){ 
			message("legend scale cannot be longer than total tree length; resetting")
			legend<-0.5*max(H)
		}
		layout(c(1,2),heights=c(0.92,0.08))
		if(outline){
			par(col="white")
			plotTree(tree,fsize=fsize[1],mar=c(0,0.1,0.1,0.1),lwd=lwd+2,offset=0.2*lwd/3+0.2/3,ftype=ftype[1])
			par(col="black")
			plotSimmap(tree,cols,pts=FALSE,lwd=lwd,fsize=fsize[1],mar=c(0,0.1,0.1,0.1),add=TRUE,ftype=ftype[1])
		} else
			plotSimmap(tree,cols,pts=FALSE,lwd=lwd,fsize=fsize[1],mar=c(0,0.1,0.1,0.1),ftype=ftype[1])
		X<-cbind(0:1000/1001,1:1001/1001)*(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label)))
		Y<-cbind(rep(0,1001),rep(0,1001))
		par(mar=c(0.1,0.1,0,0.1),xpd=NA)
		plot(NA,xlim=c(0,1),ylim=c(-0.3,0.3),xaxt="n",yaxt="n",bty="n")
		lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=lwd+2,lend=2)
		for(i in 1:1001) lines(X[i,],Y[i,],col=cols[i],lwd=lwd,lend=2)
		legf<-match(ftype[2],c("reg","b","i","bi"))
		text(x=0,y=0,round(lims[1],sig),pos=3,cex=fsize[2],font=legf)
		text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label))),y=0,round(lims[2],sig),pos=3,cex=fsize[2],font=legf)
		text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label)))/2,y=0,"trait value",pos=3,cex=fsize[2],font=legf)
		text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label)))/2,y=0,paste("length=",round(legend,3),sep=""),pos=1,cex=fsize[2],font=legf)
	} else {
		if(outline){
			par(col="white")
			plotTree(tree,cols,pts=FALSE,lwd=lwd+2,fsize=fsize[1],offset=0.2*lwd/3+0.2/3,ftype=ftype[1])
			par(col="black")
			plotSimmap(tree,cols,pts=FALSE,lwd=lwd,fsize=fsize[1],add=TRUE,ftype=ftype[1])
		} else
			plotSimmap(tree,cols,pts=FALSE,lwd=lwd,fsize=fsize[1],ftype=ftype[1])
	}
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

