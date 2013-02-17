# function plots posterior density of mapped states from stochastic mapping
# written by Liam J. Revell 2012

densityMap<-function(trees,res=100,fsize=NULL,ftype=NULL,lwd=3,check=FALSE,legend=NULL,outline=FALSE){
	tol<-1e-10
	if(class(trees)!="multiPhylo") stop("trees not 'multiPhylo' object; just use plotSimmap")
	h<-sapply(trees,function(x) max(nodeHeights(x)))
	steps<-0:res/res*max(h)
	trees<-rescaleSimmap(trees,totalDepth=max(h))
	if(check){
		X<-matrix(FALSE,length(trees),length(trees))
		for(i in 1:length(trees)) X[i,]<-sapply(trees,all.equal.phylo,current=trees[[i]])
		if(!all(X)) stop("some of the trees don't match in topology or relative branch lengths")
	}
	tree<-trees[[1]]
	H<-nodeHeights(tree)
	message("sorry - this might take a while; please be patient")
	for(i in 1:nrow(tree$edge)){
		YY<-cbind(c(H[i,1],steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))]),
			c(steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))],H[i,2]))-H[i,1]
		ZZ<-rep(0,nrow(YY))
		for(j in 1:length(trees)){
			XX<-matrix(0,length(trees[[j]]$maps[[i]]),2,dimnames=list(names(trees[[j]]$maps[[i]]),c("start","end")))
			XX[1,2]<-trees[[j]]$maps[[i]][1]
			if(length(trees[[j]]$maps[[i]])>1){
				for(k in 2:length(trees[[j]]$maps[[i]])){
					XX[k,1]<-XX[k-1,2]
					XX[k,2]<-XX[k,1]+trees[[j]]$maps[[i]][k]
				}
			}
			for(k in 1:nrow(YY)){
				lower<-which(XX[,1]<=YY[k,1]); lower<-lower[length(lower)]
				upper<-which(XX[,2]>=(YY[k,2]-tol))[1]; AA<-0
				names(lower)<-names(upper)<-NULL
				for(l in lower:upper) 
					AA<-AA+(min(XX[l,2],YY[k,2])-max(XX[l,1],YY[k,1]))/(YY[k,2]-YY[k,1])*as.numeric(rownames(XX)[l])
				ZZ[k]<-ZZ[k]+AA/length(trees)
			}
		}
		tree$maps[[i]]<-YY[,2]-YY[,1]
		names(tree$maps[[i]])<-round(ZZ*1000)
	}
	cols<-rainbow(1001,start=0.7,end=0); names(cols)<-0:1000
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
		text(x=0,y=0,"0",pos=3,cex=fsize[2],font=legf)
		text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label))),y=0,"1",pos=3,cex=fsize[2],font=legf)
		text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label)))/2,y=0,"PP(state=1)",pos=3,cex=fsize[2],font=legf)
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

