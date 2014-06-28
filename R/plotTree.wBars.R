## function to plot bars at the tips of a plotted tree
## written by Liam J. Revell 2014

plotTree.wBars<-function(tree,x,scale=1,width=NULL,type="phylogram",method="plotTree",...){
	lims<-c(-max(nodeHeights(tree))-scale*(max(x)-min(0,min(x))),max(nodeHeights(tree))+scale*(max(x)-min(0,min(x))))
	if(type=="phylogram"){
		if(method=="plotTree") capture.output(obj<-plotTree(tree,ftype="off",xlim=c(0,lims[2]),...))
		else if(method=="plotSimmap") capture.output(obj<-plotSimmap(tree,ftype="off",xlim=c(0,lims[2]),...))
	} else if(type=="fan"){
		if(method=="plotTree") capture.output(obj<-plotTree(tree,type="fan",ftype="off",xlim=lims,ylim=lims,...))
		else if(method=="plotSimmap") capture.output(obj<-plotSimmap(tree,type="fan",ftype="off",xlim=lims,ylim=lims,...))
	}
	x<-x[tree$tip.label]*scale
	if(is.null(width))
		width<-if(type=="fan") (par()$usr[4]-par()$usr[3])/(max(c(max(x)/max(nodeHeights(tree)),1))*length(tree$tip.label)) 
			else if(type=="phylogram") (par()$usr[4]-par()$usr[3])/(2*length(tree$tip.label))
	w<-width
	if(type=="phylogram"){
		if(hasArg(direction)) direction<-list(...)$direction
		else direction<-"rightwards"
		for(i in 1:length(x)){
			dx<-if(min(x)>=0) obj$xx[i] else max(obj$xx)
			dy<-obj$yy[i]
			x1<-x2<-dx
			x3<-x4<-x1+x[i]
			y1<-y4<-dy-w/2
			y2<-y3<-dy+w/2
			polygon(c(x1,x2,x3,x4)-min(0,min(x)),c(y1,y2,y3,y4),col="grey")
		}
	} else if(type=="fan"){
		if(min(x)<0) h<-max(nodeHeights(tree))
		for(i in 1:length(x)){
			theta<-atan(obj$yy[i]/obj$xx[i])
			s<-if(obj$xx[i]>0) 1 else -1
			if(min(x)>=0){
				dx<-obj$xx[i]
				dy<-obj$yy[i]
			} else {
				dx<-s*h*cos(theta)
				dy<-s*h*sin(theta)
			}
			x1<-dx+(w/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
			y1<-dy-(w/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
			x2<-dx-(w/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
			y2<-dy+(w/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
			x3<-s*x[i]*cos(theta)+x2
			y3<-s*x[i]*sin(theta)+y2
			x4<-s*x[i]*cos(theta)+x1
			y4<-s*x[i]*sin(theta)+y1
			polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col="grey")
		}
	}
	invisible(obj)
}

