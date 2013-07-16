# function plots posterior density of mapped states from stochastic mapping
# written by Liam J. Revell 2012-2013

densityMap<-function(trees,res=100,fsize=NULL,ftype=NULL,lwd=3,check=FALSE,legend=NULL,outline=FALSE,type="phylogram",plot=TRUE,...){
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	tol<-1e-10
	if(class(trees)!="multiPhylo") stop("trees not 'multiPhylo' object; just use plotSimmap")
	h<-sapply(unclass(trees),function(x) max(nodeHeights(x)))
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
	trees<-unclass(trees)
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
	x<-list(tree=tree,cols=cols); class(x)<-"densityMap"
	if(plot) plot.densityMap(x,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,type=type,mar=mar)
	invisible(x)
}

# function
# written by Liam J. Revell 2012-2013

plot.densityMap<-function(x,...){
	if(class(x)=="densityMap"){
		tree<-x$tree
		cols<-x$cols
	} else stop("x should be an object of class 'densityMap'")
	H<-nodeHeights(tree)
	# get & set optional arguments
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-FALSE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-3
	if(hasArg(leg.txt)) leg.txt<-list(...)$leg.txt
	else leg.txt<-c("0","PP(state=1)","1")
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
	if(legend){
		if(legend>max(H)){ 
			message("legend scale cannot be longer than total tree length; resetting")
			legend<-0.5*max(H)
		}
	}
	if(type=="phylogram"){
		if(legend){
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
			text(x=0,y=0,leg.txt[1],pos=3,cex=fsize[2],font=legf)
			text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label))),y=0,leg.txt[3],pos=3,cex=fsize[2],font=legf)
			text(x=(legend/max(H))*(1-fsize[1]*max(strwidth(tree$tip.label)))/2,y=0,leg.txt[2],pos=3,cex=fsize[2],font=legf)
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
	} else if(type=="fan"){
		if(outline){
			par(col="white")
			invisible(capture.output(plotTree(tree,type="fan",lwd=lwd+2,mar=mar,fsize=fsize[1],ftype=ftype[1])))
			par(col="black")
		}
		invisible(capture.output(plotSimmap(tree,cols,lwd=lwd,mar=mar,fsize=fsize[1],add=outline,ftype=ftype[1],type="fan")))
		if(legend){
			ff<-function(dd){
				if(!("."%in%dd)) dig<-0
				else dig<-length(dd)-which(dd==".")
				dig
			}
			dig<-max(sapply(strsplit(leg.txt[c(1,3)],split=""),ff))
			add.color.bar(legend,cols,title=leg.txt[2],lims<-as.numeric(leg.txt[c(1,3)]),digits=dig,prompt=FALSE,x=0.9*par()$usr[1],y=0.9*par()$usr[3],fsize=fsize[2])
		}
	}
}
